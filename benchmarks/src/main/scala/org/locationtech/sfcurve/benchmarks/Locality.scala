package org.locationtech.sfcurve.benchmarks

import java.io.{FileOutputStream, FileWriter, PrintStream}

import scala.collection.mutable.{Map=>MutableMap}

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
    val maxDPlane = Math.sqrt(360.0*360.0 + 180.0*180.0)
    def dPlane: Double = {
      val dx = a.x.degrees - b.x.degrees
      val dy = a.y.degrees - b.y.degrees
      Math.sqrt(dx * dx + dy * dy) / maxDPlane
    }

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

    def dTriangle(depth: Int): Double = {
      val t0 = TriN.getTriangle(a.x.degrees, a.y.degrees, depth)
      val t1 = TriN.getTriangle(b.x.degrees, b.y.degrees, depth)
      dBitstring(t0.bitString, t1.bitString)
    }
  }

  // utility class for testing, reporting locality
  case class Distances(name: String, curve: Any) {
    val series: Seq[String] = Seq("plane", "sphere", "index")
    val data: Map[String, MutableMap[String, Int]] = series.map(s => s -> MutableMap.empty[String, Int]).toMap

    private def acc(label: String, distance: Double): Unit = {
      require(data.contains(label))
      val submap: MutableMap[String, Int] = data(label)
      val dist = f"${distance}%1.5f"
      submap.put(dist, submap.getOrElse(dist, 0) + 1)
    }

    def accumulate(points: Points): Unit = {
      acc("plane", points.dPlane)
      acc("sphere", points.dSphere)
      val indexDistance = curve match {
        case sfc: SpaceFillingCurve2D => points.dIndex(sfc)
        case tri: Triangle => points.dTriangle(tri.depth)
      }
      acc("index", indexDistance)
    }

    def summarize(ps: PrintStream): Unit =
    {
      println(s"\n$name distances...")
      series.foreach(label => {
        println(s"  '$label':")
        val dists: MutableMap[String, Int] = data(label)
        val nTransitions = dists.values.sum
        println(s"    Number of cell-to-cell transitions:  $nTransitions")
        val totalDist = dists.map {
          case (k, v) => k.toDouble * v
        }.sum
        println(s"    Mean inter-cell distance:  ${totalDist / nTransitions.toDouble}")
      })
    }

  }

  def generateLongitude: Degrees = Degrees(0.5 + Math.random() * 359.0 - 180.0)

  def generateLatitude: Degrees = Degrees(0.5 + Math.random() * 179.0 - 90.0)

  def generatePoint: Point = Point(generateLongitude, generateLatitude)

  def generatePoints: Points = Points(generatePoint, generatePoint)

  case class CellIterator(name: String, curve: {def cardinality: Long}) {
    val cardinality: Long = curve.cardinality
    val distances: Distances = Distances(name, curve)

    var lastState: Any = _
    var state: Any = _
    var counter: Int = 0

    state = curve match {
      case sfc: SpaceFillingCurve2D => 0
      case tri: Triangle => tri
    }
    var lastPoint: Point = getPoint
    counter = 1
    lastState = state
    state = curve match {
      case sfc: SpaceFillingCurve2D => 1
      case tri: Triangle => tri.next.orNull
    }
    var point: Point = getPoint

    private def getPoint: Point = curve match {
      case sfc: SpaceFillingCurve2D =>
        val p = sfc.toPoint(counter)
        Point(Degrees(p._1), Degrees(p._2))
      case tri: Triangle =>
        val result = Point(Degrees(state.asInstanceOf[Triangle].xMid), Degrees(state.asInstanceOf[Triangle].yMid))
        result
    }

    def hasNext: Boolean = counter < cardinality

    def next(): Points = {
      val result = Points(lastPoint, point)
      //println(s"next:  result = $result")
      lastPoint = point
      lastState = state
      counter = counter + 1
      if (hasNext) {
        state = curve match {
          case sfc: SpaceFillingCurve2D => counter
          case tri: Triangle => state.asInstanceOf[Triangle].next.orNull
        }
        point = getPoint
      }
      result
    }

    def exhaust(ps: PrintStream): Unit = {
      while (hasNext) {
        distances.accumulate(next())
      }
      distances.summarize(ps)
    }
  }

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
    val TriangleDepth = bitsPrecision.toInt / 3

    CellIterator("Z2", z2).exhaust(ps)
    CellIterator("H2", h2).exhaust(ps)
    CellIterator("triangle", TriN.createLowestIndex(TriangleDepth)).exhaust(ps)
  } catch { case t: Throwable =>
    ps.close()
  }

}
