package org.locationtech.sfcurve.benchmarks

import java.io.{FileOutputStream, FileWriter, PrintStream}

import org.locationtech.sfcurve.hilbert.HilbertCurve2D
import org.locationtech.sfcurve.zorder.{TriN, Triangle, ZCurve2D}
import org.locationtech.sfcurve.{Dimensions, SpaceFillingCurve2D}

object Locality extends App {
  def toBinaryString(index: Long, bits: Int): String =
    index.toBinaryString.reverse.padTo(bits, "0").reverse.mkString("")

  def dBitstring(idxA: String, idxB: String): Double = {
    val result: Double = idxA.zip(idxB).zipWithIndex.map {
      case ((bitA, bitB), pos) if bitA != bitB =>
        //println(s"  bits mismatch:  $bitA, $bitB")
        Math.pow(0.5, pos + 1)
      case ((bitA, bitB), pos) if bitA == bitB =>
        //println(s"  bits MATCH:  $bitA, $bitB")
        0.0
    }.sum

    //println(s"dBitString($idxA, $idxB) = $result")

    result
  }

  case class Degrees(degrees: Double) {
    def radians: Double = Math.toRadians(degrees)
  }

  case class Point(x: Degrees, y: Degrees) {
    def wkt: String = "POINT(" + x.degrees + " " + y.degrees + ")"
  }

  case class Points(a: Point, b: Point) {
    def dSphere: Double = {
      val x0 = a.x.radians
      val y0 = a.y.radians
      val x1 = b.x.radians
      val y1 = b.y.radians
      Math.acos(Math.sin(y0) * Math.sin(y1) + Math.cos(y0) * Math.cos(y1) * Math.cos(Math.abs(x1 - x0))) / Math.PI
    }

    implicit def deg2dec(d: Degrees): Double = d.degrees

    implicit def dec2deg(d: Double): Degrees = Degrees(d)

    def dIndex(curve: SpaceFillingCurve2D): Double = {
      val bitsPrecision: Int = Dimensions.bitsFromCardinality(curve.cardinality).toInt
      val idxA: String = toBinaryString(curve.toIndex(a.x, a.y), bitsPrecision)
      val idxB: String = toBinaryString(curve.toIndex(b.x, b.y), bitsPrecision)
      dBitstring(idxA, idxB)
    }

    def dTriangle(bits: Int): Double = {
      val t0 = TriN.getTriangle(a.x.degrees, a.y.degrees, bits)
      val t1 = TriN.getTriangle(b.x.degrees, b.y.degrees, bits)
      dBitstring(t0.bitString, t1.bitString)
    }
  }

  def generateLongitude: Degrees = Degrees(0.5 + Math.random() * 359.0 - 180.0)

  def generateLatitude: Degrees = Degrees(0.5 + Math.random() * 179.0 - 90.0)

  def generatePoint: Point = Point(generateLongitude, generateLatitude)

  def generatePoints: Points = Points(generatePoint, generatePoint)

  // set up
  val bitsPrecision: Long = 6
  val cardinality = 1L << bitsPrecision
  val cardPerDim = 1L << (bitsPrecision >> 1L)
  println(s"Bits precision $bitsPrecision, cardinality $cardinality")
  val z2 = new ZCurve2D(cardPerDim.toInt)
  println(s"Z2 cardinality ${z2.cardinality}")
  val h2 = new HilbertCurve2D(cardPerDim.toInt)
  println(s"H2 cardinality ${h2.cardinality}")

  // unit testing
  val point = Point(Degrees(-78.0), Degrees(38.0))
  val points = Points(point, point)
  val dIndexZ2 = points.dIndex(z2)
  println(s"dIndex(z2):  $dIndexZ2")
  assert(dIndexZ2 == 0.0)
  val dIndexH2 = points.dIndex(h2)
  println(s"dIndex(h2):  $dIndexH2")
  assert(dIndexH2 == 0.0)
  val dSphere = points.dSphere
  println(s"dSphere(z2):  $dSphere")
  assert(dSphere >= 0.0 && dSphere <= 0.00001)
  assert(dBitstring("100", "100") == 0.0)
  assert(dBitstring("110", "100") == 0.25)
  assert(dBitstring("100", "101") == 0.125)

  // application
  // =CORREL($C$2:$C$100001,D$2:D$100001)
  // Z2	H2	tri
  //0.716675856	0.250055058	0.743290119
  //0.714665639	0.240273265	0.741659422
  //0.712142341	0.242960041	0.740031498
  //0.719532321	0.247919821	0.745290572
  //0.713298924	0.237534539	0.740847436
  //0.714175341	0.242465329	0.741506736
  //0.715370096	0.248664902	0.742935838
  //0.717145868	0.242575736	0.743071597
  //0.715021582	0.239930681	0.741891317
  //0.715148242	0.238397849	0.740885343


  var ps: PrintStream = null
  try {
    /*
    ps = new PrintStream(new FileOutputStream("test-distance.csv"))
    val curves = Seq(z2, h2)
    //val ps: PrintStream = System.out
    ps.println("point_a_wkt,point_b_wkt,dist_sphere,dist_z2,dist_h2,dist_tri")
    for (i <- 1 to 100000) {
      val points = generatePoints
      ps.print(s"""\"${points.a.wkt}\",\"${points.b.wkt}\"""")
      ps.print(s",${points.dSphere}")
      for (sfc <- curves) {
        //ps.print(s"\nindex A:  ${toBinaryString(sfc.toIndex(points.a.x.degrees, points.a.y.degrees), bitsPrecision.toInt)}")
        //ps.print(s"\nindex B:  ${toBinaryString(sfc.toIndex(points.b.x.degrees, points.b.y.degrees), bitsPrecision.toInt)}")
        ps.print(s",${points.dIndex(sfc)}")
      }
      ps.print(s",${points.dTriangle(bitsPrecision.toInt / 3)}")
      ps.println()
    }
    */

    var z2Idx = 0L
    val z2Dists = collection.mutable.Map.empty[String, Int]
    while (z2Idx < z2.cardinality) {
      if (z2Idx > 0) {
        val xy0 = z2.toPoint(z2Idx - 1)
        val a = Point(Degrees(xy0._1), Degrees(xy0._2))
        val xy1 = z2.toPoint(z2Idx)
        val b = Point(Degrees(xy1._1), Degrees(xy1._2))
        val points = Points(a, b)
        val dist = f"${points.dSphere}%1.5f"
        z2Dists.put(dist, z2Dists.getOrElse(dist, 0) + 1)
      }
      z2Idx = z2Idx + 1
    }
    {
      println("\nZ2 distances...")
      println(s"  Total inter-Z2 distances:  ${z2Dists.values.sum}")
      val totalDist = z2Dists.map {
        case (k, v) => k.toDouble * v
      }.sum
      println(s"  Mean inter-Z2 distance:  ${totalDist / z2Dists.values.sum}")
    }

    var h2Idx = 0L
    val h2Dists = collection.mutable.Map.empty[String, Int]
    while (h2Idx < h2.cardinality) {
      if (h2Idx > 0) {
        val xy0 = h2.toPoint(z2Idx - 1)
        val a = Point(Degrees(xy0._1), Degrees(xy0._2))
        val xy1 = h2.toPoint(z2Idx)
        val b = Point(Degrees(xy1._1), Degrees(xy1._2))
        val points = Points(a, b)
        val dist = f"${points.dSphere}%1.5f"
        h2Dists.put(dist, h2Dists.getOrElse(dist, 0) + 1)
      }
      h2Idx = h2Idx + 1
    }
    {
      println("\nH2 distances...")
      println(s"  Total inter-H2 distances:  ${h2Dists.values.sum}")
      val totalDist = h2Dists.map {
        case (k, v) => k.toDouble * v
      }.sum
      println(s"  Mean inter-H2 distance:  ${totalDist / h2Dists.values.sum}")
    }

    val depth: Int = bitsPrecision.toInt / 3
    val triCardinality = 8 * (depth - 1) * 4
    var triOpt: Option[Triangle] = Option(TriN.invIndex(List.fill(depth)(TriN.TransCenter).foldLeft(0L)((acc, t) => (acc << 3L) | t), depth))
    var lastTriOpt: Option[Triangle] = None
    val triDists = collection.mutable.Map.empty[String, Int]
    while (triOpt.isDefined) {
      println(s"Triangle ${triOpt.get.bitString}, last $lastTriOpt...")
      lastTriOpt.foreach(lastTri => {
        val point0 = Point(Degrees(lastTri.xMid), Degrees(lastTri.yMid))
        val point1 = Point(Degrees(triOpt.get.xMid), Degrees(triOpt.get.yMid))
        val points = Points(point0, point1)
        val dist = f"${points.dSphere}%1.5f"
        triDists.put(dist, triDists.getOrElse(dist, 0) + 1)
      })
      lastTriOpt = triOpt
      triOpt = triOpt.get.next
      println(s"got next $triOpt")
    }
    {
      println("\nTriangle distances...")
      println(s"  Total inter-triangle distances:  ${triDists.values.sum}")
      val totalDist = triDists.map {
        case (k, v) => k.toDouble * v
      }.sum
      println(s"  Mean inter-triangle distance:  ${totalDist / triDists.values.sum}")
      println(s"  Triangle cardinality:  $triCardinality")
    }

  } catch { case t: Throwable =>
    ps.close()
  }

}
