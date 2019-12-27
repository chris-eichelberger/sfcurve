package org.locationtech.sfcurve.benchmarks

import java.io.{FileOutputStream, FileWriter, PrintStream}

import scala.collection.mutable.{Map => MutableMap}
import org.locationtech.sfcurve.hilbert.HilbertCurve2D
import org.locationtech.sfcurve.zorder.{TriN, Triangle, ZCurve2D}
import org.locationtech.sfcurve.{Dimensions, IndexRange, RangeComputeHints, SpaceFillingCurve2D}

object Locality extends App {
  type Curve = {
    def cardinality: Long;
    def toIndex(x: Double, y: Double): Long;
    def toPoint(index: Long): (Double, Double);
    def toRanges(xmin: Double, ymin: Double, xmax: Double, ymax: Double, hints: Option[RangeComputeHints]): Seq[IndexRange]
  }

  def toBinaryString(index: Long, bits: Int): String =
    index.toBinaryString.reverse.padTo(bits, "0").reverse.mkString("")

  def dBitstring(idxA: String, idxB: String, verbose: Boolean = false): Double = {
    val result: Double = idxA.zip(idxB).zipWithIndex.map {
      case ((bitA, bitB), pos) if bitA != bitB =>
        //println(s"  bits mismatch:  $bitA, $bitB")
        Math.pow(0.5, pos + 1)
      case ((bitA, bitB), pos) if bitA == bitB =>
        //println(s"  bits MATCH:  $bitA, $bitB")
        0.0
    }.sum

    if (verbose) {
      println(s"dBitString($idxA, $idxB) = $result")
    }

    result
  }

  case class Degrees(degrees: Double) {
    def radians: Double = Math.toRadians(degrees)
  }

  case class Point(x: Degrees, y: Degrees) {
    def wkt: String = "POINT(" + x.degrees + " " + y.degrees + ")"
  }

  case class Points(a: Point, b: Point, verbose: Boolean = false) {
    val maxDPlane = Math.sqrt(360.0*360.0 + 180.0*180.0)
    def dPlane(verbose: Boolean = false): Double = {
      val dx = a.x.degrees - b.x.degrees
      val dy = a.y.degrees - b.y.degrees
      val result = Math.sqrt(dx * dx + dy * dy) / maxDPlane
      if (verbose) {
        println(s"dPlane:  ${a.wkt} to ${b.wkt} = $result")
      }
      result
    }

    def dSphere(verbose: Boolean = false): Double = {
      val x0 = a.x.radians
      val y0 = a.y.radians
      val x1 = b.x.radians
      val y1 = b.y.radians
      val result = Math.acos(Math.sin(y0) * Math.sin(y1) + Math.cos(y0) * Math.cos(y1) * Math.cos(Math.abs(x1 - x0))) / Math.PI
      if (verbose) {
        println(s"dSphere:  ${a.wkt} to ${b.wkt} = $result")
      }
      result
    }

    implicit def deg2dec(d: Degrees): Double = d.degrees

    implicit def dec2deg(d: Double): Degrees = Degrees(d)

    def dIndex(curve: SpaceFillingCurve2D, verbose: Boolean = false): Double = {
      val bitsPrecision: Int = Dimensions.bitsFromCardinality(curve.cardinality).toInt
      val idxA: String = toBinaryString(curve.toIndex(a.x, a.y), bitsPrecision)
      val idxB: String = toBinaryString(curve.toIndex(b.x, b.y), bitsPrecision)
      dBitstring(idxA, idxB, verbose = verbose)
    }

    def dTriangle(depth: Int, verbose: Boolean = false): Double = {
      val t0 = TriN.getTriangle(a.x.degrees, a.y.degrees, depth)
      val t1 = TriN.getTriangle(b.x.degrees, b.y.degrees, depth)
      dBitstring(t0.bitString, t1.bitString, verbose = verbose)
    }
  }

  // utility class for testing, reporting locality
  case class Distances(name: String, curve: {def toIndex(x: Double, y: Double): Long; def toPoint(i: Long): (Double, Double)}) {
    val series: Seq[String] = Seq("plane", "sphere", "index")
    val data: Map[String, MutableMap[String, Int]] = series.map(s => s -> MutableMap.empty[String, Int]).toMap

    private def acc(label: String, distance: Double): Unit = {
      require(data.contains(label))
      val submap: MutableMap[String, Int] = data(label)
      val dist = f"${distance}%1.5f"
      submap.put(dist, submap.getOrElse(dist, 0) + 1)
    }

    def cellPoints(points: Points): Points = {
      val ai2: Long = curve.toIndex(points.a.x.degrees, points.a.y.degrees)
      val a2: Point = Point(Degrees(curve.toPoint(ai2)._1), Degrees(curve.toPoint(ai2)._2))
      val bi2: Long = curve.toIndex(points.b.x.degrees, points.b.y.degrees)
      val b2: Point = Point(Degrees(curve.toPoint(bi2)._1), Degrees(curve.toPoint(bi2)._2))
      Points(a2, b2)
    }

    def accumulate(rawPoints: Points, verbose: Boolean = false): Unit = {
      val points: Points = cellPoints(rawPoints)
      acc("plane", points.dPlane(verbose = verbose))
      acc("sphere", points.dSphere(verbose = verbose))
      val indexDistance = curve match {
        case sfc: SpaceFillingCurve2D => points.dIndex(sfc, verbose = verbose)
        case tri: Triangle => points.dTriangle(tri.depth, verbose = verbose)
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

  trait Sampler extends Iterator[Points] {
    def name: String
    def curve: Curve
  }

  case class CellIterator(name: String, curve: Curve) extends Sampler {
    val cardinality: Long = curve.cardinality

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
        require(state != null, s"'state' was <NULL> unexpectedly for counter $counter (< cardinality $cardinality)")
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
  }

  case class FixedSample(name: String, curve: Curve, input: Seq[Points]) extends Sampler {
    val itr: Iterator[Points] = input.iterator
    def hasNext: Boolean = itr.hasNext
    def next(): Points = itr.next()
  }

  case class RandomSample(name: String, curve: Curve, numPoints: Long) extends Sampler {
    var counter: Long = 0
    def hasNext: Boolean = counter < numPoints
    def next(): Points = { counter = counter + 1; generatePoints }
  }

  case class Summary(sampler: Sampler, verbose: Boolean = false) {
    val distances: Distances = Distances(sampler.name, sampler.curve)

    def exhaust(ps: PrintStream): Unit = {
      while (sampler.hasNext) {
        distances.accumulate(sampler.next(), verbose = verbose)
      }
      distances.summarize(ps)
    }
  }

  // set up
  val bitsPrecision: Long = 54
  require((bitsPrecision % 2) == 0, "bitsPrecision must be divisible by 2 for Z2, H2")
  require((bitsPrecision % 3) == 0, "bitsPrecision must be divisible by 3 for T2")
  val cardinality = 1L << bitsPrecision
  val cardPerDim = 1L << (bitsPrecision >> 1L)
  println(s"Bits precision $bitsPrecision, cardinality $cardinality")
  val z2 = new ZCurve2D(cardPerDim.toInt)
  println(s"Z2 cardinality ${z2.cardinality}")
  val h2 = new HilbertCurve2D(cardPerDim.toInt)
  println(s"H2 cardinality ${h2.cardinality}")
  val TriangleDepth = bitsPrecision.toInt / 3
  val t2 = TriN.createLowestIndex(TriangleDepth)
  println(s"Triangle depth:  ${TriangleDepth}")
  println(s"Triangle cardinality:  ${t2.cardinality}")

  // unit testing
  val point = Point(Degrees(-78.0), Degrees(38.0))
  val points = Points(point, point)
  val dIndexZ2 = points.dIndex(z2)
  println(s"dIndex(z2):  $dIndexZ2")
  assert(dIndexZ2 == 0.0)
  val dIndexH2 = points.dIndex(h2)
  println(s"dIndex(h2):  $dIndexH2")
  assert(dIndexH2 == 0.0)
  val dSphere = points.dSphere()
  println(s"dSphere(z2):  $dSphere")
  assert(dSphere >= 0.0 && dSphere <= 0.00001)
  assert(dBitstring("100", "100") == 0.0)
  assert(dBitstring("110", "100") == 0.25)
  assert(dBitstring("100", "101") == 0.125)
  assert(TriN.createLowestIndex(1).cardinality == 8)
  assert(TriN.createLowestIndex(2).cardinality == 32)
  assert(TriN.createLowestIndex(3).cardinality == 128)
  assert(TriN.createLowestIndex(4).cardinality == 512)

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

  val exhaustive = Seq[Sampler](
    CellIterator("Z2", z2),
    CellIterator("H2", h2),
    CellIterator("T2", t2)
  )

  val numRandomPoints: Long = 10000
  val random = Seq[Sampler](
    RandomSample("Z2", z2, numRandomPoints),
    RandomSample("H2", h2, numRandomPoints),
    RandomSample("T2", t2, numRandomPoints)
  )

  val CharlottesvillePoints = Seq[Points](
    Points(
      Point(Degrees(-78.495150), Degrees(38.075776)),  // CCRi
      Point(Degrees(-78.688256), Degrees(38.054444))  // Eichelberger
    )
  )
  val cville = Seq[Sampler](
    FixedSample("Z2", z2, CharlottesvillePoints),
    FixedSample("H2", h2, CharlottesvillePoints),
    FixedSample("T2", t2, CharlottesvillePoints)
  )

  val samples: Seq[Sampler] = cville

  for (sample <- samples) {
    Summary(sample, verbose = true).exhaust(ps)
  }

//  try {
//  } catch { case t: Throwable =>
//    ps.close()
//  }

}
