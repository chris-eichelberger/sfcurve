package org.locationtech.sfcurve.zorder

import java.io.{FileWriter, PrintWriter}

import org.locationtech.sfcurve.Dimensions.{Cell, Extent, Latitude, Longitude}
import org.locationtech.sfcurve.{InMemoryRangeConsolidator, IndexRange, RangeComputeHints, SpaceFillingCurve2D}
import org.locationtech.sfcurve.Utilities.CartesianProductIterable

/**
  * Definition of triangle shapes / orientations:
  * Apex on top...
  *   clockwise...
  *     0 - A|BC
  *     1 - B|CA
  *     2 - C|AB
  *   counter clockwise...
  *     3 - A|CB
  *     4 - B|AC
  *     5 - C|BA
  * Apex on bottom
  *   clockwise...
  *     6 - AB|C
  *     7 - BC|A
  *     8 - CA|B
  *   counter clockwise...
  *     9 - AC|B
  *     A - BA|C
  *     B - CB|A
  */

import TriN._

case class Rectangle(x: Extent[Double], y: Extent[Double]) {
  def isPoint: Boolean = x.isPoint && y.isPoint
}

// although this will often have a degenerate Apex Geo X (when the shape is a proper triangle),
// it can also describe a real rectangle (as for an octahedron face) as well as a trapezoid
case class TriBounds(octXApex: Double, octXBase: Extent[Double], Y: Extent[Double], isApexUp: Boolean) {
  val octX0: Double = octXBase.min
  val octX1: Double = octXBase.max
  val octXMid: Double = octXBase.q2
  val octX14: Double = octXBase.q1
  val octX34: Double = octXBase.q3

  val (yApex: Double, yBase: Double) = if (isApexUp) (Y.max, Y.min) else (Y.min, Y.max)

  // TODO fix this so that it doesn't return two copies of the same point at the singularity
  val geoXApex: Extent[Double] = Extent(TriN.octToGeoX(yApex)(octXApex - ApexBias), TriN.octToGeoX(yApex)(octXApex + ApexBias))
  val geoXBase: Extent[Double] = octXBase.apply(TriN.octToGeoX(yBase))

  val geoX0: Double = geoXBase.min
  val geoX1: Double = geoXBase.max
  val geoXMid: Double = geoXBase.q2
  val geoX14: Double = geoXBase.q1
  val geoX34: Double = geoXBase.q3

  val y0: Double = Y.min
  val y1: Double = Y.max
  val yMid: Double = Y.q2

  def geoCenter: (Double, Double) = (geoXMid, yMid)

  def octCenter: (Double, Double) = (octXMid, yMid)

  def geoXAtMidY: Extent[Double] = Extent(
    0.5 * (geoXApex.min + geoXBase.min),
    0.5 * (geoXApex.max + geoXBase.max)
  )

  def octContained(octXMin: Double, ymin: Double, octXMax: Double, ymax: Double): Boolean = {
    octXMin <= octXBase.min && octXMax >= octXBase.max && ymin <= Y.min && ymax >= Y.max
  }

  def octContainsPoint(octX: Double, y: Double): Boolean = {
    // simple rejection
    if (!octXBase.contains(octX)) return false
    if (!Y.contains(y)) return false

    val rawPY: Double = (y - y0) / (y1 - y0)
    val pY: Double = if (isApexUp) 1.0 - rawPY else rawPY
    val halfWidth: Double = octXBase.span * 0.5 * pY
    val dist: Double = Math.abs(octX - octXBase.q2)

    // TODO remove after debugging
//    println(f"octContainsPoint($octX%1.3f, $y%1.3f):  |$octX%1.3f - ${octXBase.q2}%1.3f| = $dist%1.3f <= $halfWidth?  ${dist <= halfWidth}%s")
//    println(f"  isApexUp?  $isApexUp%s")
//    println(f"  rawPY ($y%1.3f - $y0%1.3f) / ($y1%1.3f - $y0%1.3f) = ${y - y0}%1.3f / ${y1 - y0}%1.3f = $rawPY%1.3f")
//    println(f"  pY $pY%1.3f")

    dist <= halfWidth
  }

  def octOverlaps(rectangle: Rectangle): Boolean = {
    // TODO:  remove after debugging
//    val r: String = f"R(X ${rectangle.x.min}%1.3f, ${rectangle.x.max}%1.3f Y ${rectangle.y.min}%1.3f, ${rectangle.y.max}%1.3f)"
//    val t: String = f"T(X ${octXBase.min}%1.3f, ${octXBase.max}%1.3f Y ${Y.min}%1.3f ${Y.max}%1.3f apex $octXApex%1.3f, up $isApexUp%s)"

    // point elimination
    if (rectangle.isPoint) return octContainsPoint(rectangle.x.min, rectangle.y.min)

    // simple vertical elimination
    if (rectangle.y.max < Y.min) {
      // TODO:  remove after debugging
//      println(f"    TriBounds.overlaps $r%s, $t%s :  r.y.max ${rectangle.y.max}%1.3f < Y.min ${Y.min}%1.3f")
      return false
    }
    if (rectangle.y.min > Y.max) {
      // TODO:  remove after debugging
//      println(f"    TriBounds.overlaps $r%s, $t%s :  r.y.min ${rectangle.y.min}%1.3f > Y.max ${Y.max}%1.3f")
      return false
    }

    // simple horizontal elimination
    if (rectangle.x.max < octXBase.min) {
      // TODO:  remove after debugging
//      println(f"    TriBounds.overlaps $r%s, $t%s :  r.x.max ${rectangle.x.max}%1.3f < octXBase.min ${octXBase.min}%1.3f")
      return false
    }
    if (rectangle.x.min > octXBase.max) {
      // TODO:  remove after debugging
//      println(f"    TriBounds.overlaps $r%s, $t%s :  r.x.min ${rectangle.x.min}%1.3f > octXBase.max ${octXBase.max}%1.3f")
      return false
    }

    // simple acceptance:  rectangle contains triangle corners
    if (rectangle.x.contains(octXApex) && rectangle.y.overlaps(Y)) {
      // TODO:  remove after debugging
//      println(f"    Simple RcT acceptance of apex")
      return true
    }
    if (rectangle.y.contains(yBase) && rectangle.x.contains(octX0)) {
      // TODO:  remove after debugging
//      println(f"  Simple RcT acceptance of minX corner")
      return true
    }
    if (rectangle.y.contains(yBase) && rectangle.x.contains(octX1)) {
      // TODO:  remove after debugging
//      println(f"  Simple RcT acceptance of maxX corner")
      return true
    }

    // simple acceptance:  triangle contains rectangle points
    for (rX <- Seq(rectangle.x.q0, rectangle.x.q2, rectangle.x.q4); rY <- Seq(rectangle.y.q0, rectangle.y.q2, rectangle.y.q4)) {
      if (octContainsPoint(rX, rY)) {
        // TODO:  remove after debugging
//        println(f"    Simple TcR")
        return true
      }
    }

    // shoulder cases, both left and right
    val slope: Double = if (rectangle.x.max <= octXBase.q2 ^ isApexUp) MNeg else MPos
    val b: Double = if (isApexUp) Y.min else Y.max
    val x0: Double = if (rectangle.x.max <= octXBase.q2) octXBase.min else octXBase.max
    val (xProbe: Double, yProbe: Double) =
      if (rectangle.x.contains(octXBase.q2)) {
        // TODO remove after debugging!
//        println(f"    xProbe:  Using octXBase.q2 ${octXBase.q2}%1.3f contained within R")
//        println(f"      yProbe:  computed ${min(y1, rectangle.y.max, max(y0, rectangle.y.min))}")
        (octXBase.q2, min(y1, rectangle.y.max, max(y0, rectangle.y.min)))
      } else {
        if (octXBase.q2 <= rectangle.x.min) {
          // TODO remove after debugging!
//          println(f"    xProbe:  Right of oct mid (${octXBase.q2}), using r.x.min ${rectangle.x.min}%1.3f")
          if (rectangle.y.contains(y0)) {
//            println(f"      yProbe:  Using y0")
            (rectangle.x.min, y0)
          } else {
            // TODO remove after debugging!
//            println(f"    rX ${rectangle.x.toString}%s, tX ${octXBase.toString}%s")
//            println(f"    rY ${rectangle.y.toString}%s, tY ${Y.toString}%s")
            return false
            //throw new Exception("how do you get here #1?")
          }
        } else {
          // TODO remove after debugging!
//          println(f"    xProbe:  Left of oct mid (${octXBase.q2}), using r.x.max ${rectangle.x.max}%1.3f")
          if (rectangle.y.contains(y0)) {
//            println(f"      yProbe:  Using y0")
            (rectangle.x.max, y0)
          } else {
//            println(f"    rX ${rectangle.x.toString}%s, tX ${octXBase.toString}%s")
//            println(f"    rY ${rectangle.y.toString}%s, tY ${Y.toString}%s")
//            println(s"    t.contains(${rectangle.x.max}, ${rectangle.y.min}) ${octContainsPoint(rectangle.x.max, rectangle.y.min)}")
            return false
            //throw new Exception("how do you get here #2?")
          }
        }
      }
    require(xProbe >= octXBase.min && xProbe <= octXBase.max, "Invalid xProbe")
    require(yProbe >= Y.min && yProbe <= Y.max)
    val yEq = slope * (xProbe - x0) + b
    val inside = if (isApexUp) yProbe <= yEq else yProbe >= yEq

    // TODO:  remove after debugging
//    println(f"    TriBounds.overlaps($r%s, $t%s:  yProbe $yProbe%1.3f ${if(isApexUp) "<=" else ">="} $slope%1.3f * ($xProbe%1.3f - $x0%1.3f) + $b%1.3f = $yEq%1.3f, inside $inside%s")
//    println(f"      Triangle X-mid ${octXBase.q2}%1.3f")

    inside
  }

  // assumes OCTAHEDRAL coordinates
  def octOverlaps(xmin: Double, ymin: Double, xmax: Double, ymax: Double): Boolean = {
    octOverlaps(Rectangle(Extent(xmin, xmax), Extent(ymin, ymax)))
  }

  def containsOctX(x: Double): Boolean = x >= octX0 && x <= octX1

  def containsY(y: Double): Boolean = y >= y0 && y <= y1

  def areaSquareDegrees: Double = 0.5 * Math.abs(geoX1 - geoX0) * Math.abs(y1 - y0)

  def areaSquareMeters: Double = {
    val dxDeg = Math.abs(geoX1 - geoX0)
    val dyDeg = Math.abs(y1 - y0)
    val yDeg = 0.5 * (y1 + y0)
    val dyM = dyDeg * 111.32 * 1000.0
    val dxM = dxDeg * Math.cos(Math.toRadians(yDeg)) * 111.32 * 1000.0
    0.5 * dyM * dxM
  }

  def octWkt(orientation: Int): String = orientation match {
    case A_BC | A_CB | B_AC | B_CA | C_AB | C_BA =>
      // apex on top
      s"POLYGON(($octXApex ${Y.max}, ${octXBase.max} ${Y.min}, ${octXBase.min} ${Y.min}, $octXApex ${Y.max}))"
    case AB_C | AC_B | BA_C | BC_A | CA_B | CB_A =>
      // apex on bottom
      s"POLYGON((${octXBase.min} ${Y.max}, ${octXBase.max} ${Y.max}, $octXApex ${Y.min}, ${octXBase.min} ${Y.max}))"
    case _ =>
      throw new Exception(s"Invalid orientation $orientation")
  }

  def geoWkt(orientation: Int): String = orientation match {
    case A_BC | A_CB | B_AC | B_CA | C_AB | C_BA =>
      // apex on top
      s"POLYGON((${geoXApex.min} ${Y.max}, ${geoXApex.max} ${Y.max}, ${geoXBase.max} ${Y.min}, ${geoXBase.min} ${Y.min}, ${geoXApex.min} ${Y.max}))"
    case AB_C | AC_B | BA_C | BC_A | CA_B | CB_A =>
      // apex on bottom
      s"POLYGON((${geoXBase.min} ${Y.max}, ${geoXBase.max} ${Y.max}, ${geoXApex.max} ${Y.min}, ${geoXApex.min} ${Y.min}, ${geoXBase.min} ${Y.max}))"
    case _ =>
      throw new Exception(s"Invalid orientation $orientation")
  }
}

/**
  *
  * @param index the whole number uniquely representing this triangle-cell index
  * @param orientation which Gray-encoding orientation to use for this cell
  * @param bounds the geographic and octahedral bounds of this triangle
  * @param depth how many levels deep this cell is embedded
  */
case class Triangle(index: Long, orientation: Int, bounds: TriBounds, depth: Int) extends SpaceFillingCurve2D(depth * 3) with InMemoryRangeConsolidator {
  implicit def octX0: Double = bounds.octX0
  implicit def octX1: Double = bounds.octX1
  implicit def octXMid: Double = bounds.octXMid
  implicit def octX14: Double = bounds.octX14
  implicit def octX34: Double = bounds.octX34

  val (y0: Double, y1: Double, mLeftInv: Double, mRightInv: Double) = if (isApexUp(orientation)) {
    (bounds.Y.min, bounds.Y.max, MNegInv, MPosInv)
  } else {
    (bounds.Y.max, bounds.Y.min, MPosInv, MNegInv)
  }
  val yMid: Double = 0.5 * (y0 + y1)

  override val cardinality: Long = 8L * Math.pow(4L, depth.toLong - 1L).toLong

  val name: String = s"Triangle${depth}"

  def bitString: String = indexBinaryString(index, depth)

  // degenerate for the Curve contract
  // assumes that you are given GEOGRAPHIC coordinates!
  override def toIndex(geoX: Double, y: Double): Long = {
    TriN.index(geoX, y, depth)
  }

  // degenerate for the Curve contract
  override def toPoint(index: Long): (Double, Double) = {
    TriN.invIndex(index, depth).bounds.geoCenter
  }

  def expandedIndex(cIdx: Long): Long = TriN.expandedIndex(TriN.longToBinaryString(cIdx, 3 + 2 * (depth -1)))

  def compactIndex: Long = TriN.compactIndex(bitString)

  // every triangle has four smaller sub-triangle children
  def childTriangles: Seq[Triangle] = Transitions.toSeq.map(child)

  def getRangeIndex: Long = if (TriN.UseCompactIndex) compactIndex else index

  def getChildIndexRange(maxDepth: Int): IndexRange = {
    require(maxDepth >= depth)

    val firstDescendant: Triangle = TriN.invIndex(bitString + "000" * (maxDepth - depth))
    val lastDescendant: Triangle = TriN.invIndex(bitString + "100" * (maxDepth - depth))
    require(firstDescendant.index < lastDescendant.index, s"Child index-range values were returned out of order:  $bitString -> ${firstDescendant.bitString} to ${lastDescendant.bitString}")

    IndexRange(firstDescendant.getRangeIndex, lastDescendant.getRangeIndex, contained = true)
  }

  // called by TriN and itself
  // assumes that you are given OCTAHEDRAL coordinates!
  def getRangesRecursively(octXMin: Double, ymin: Double, octXMax: Double, ymax: Double, maxDepth: Int): Seq[IndexRange] = {
    // check stop conditions
    if (depth == maxDepth) {
      // if you've gotten here, then your (single) index is the only thing to return
      return Seq(IndexRange(getRangeIndex, getRangeIndex, contained = true))
    }
    if (depth > maxDepth) throw new Exception("Should not be able to recurse this far")

    // if you're entirely inside the query area, just return everything (down to the maximum depth)
    if (bounds.octContained(octXMin, ymin, octXMax, ymax))
      return Seq(getChildIndexRange(maxDepth))

    // find your children that intersect the query area
    val subs = childTriangles.filter(child => child.bounds.octOverlaps(octXMin, ymin, octXMax, ymax))
    if (subs.isEmpty) {
      throw new Exception(s"Parent\n  $this\noverlaps oct bounds\n  ($octXMin, $ymin, $octXMax, $ymax)\nbut no child does:${childTriangles.map(_.toString).mkString("\n  ", "\n  ", "")}")
    }

    // recurse, and combine the results
    subs.tail.foldLeft(subs.head.getRangesRecursively(octXMin, ymin, octXMax, ymax, maxDepth))((acc, t) => {
      t.consolidateRanges((acc ++ t.getRangesRecursively(octXMin, ymin, octXMax, ymax, maxDepth)).iterator).toSeq
    })
  }

  def parentOpt: Option[Triangle] = {
    if (depth == 1) return None
    Option(TriN.invIndex(index >> 3L, depth - 1))
  }

  // degenerate for the Curve contract
  // assumes that you are given GEOGRAPHIC coordinates!
  override def toRanges(geoXMin: Double, ymin: Double, geoXMax: Double, ymax: Double, hints: Option[RangeComputeHints] = None): Seq[IndexRange] = {
    require(geoXMin <= geoXMax)
    require(ymin <= ymax)

    // defer to the companion object to do this planning top-down
    TriN.getRanges(geoXMin, ymin, geoXMax, ymax, depth)
  }

  def fold(subordinates: Seq[Long]): Long = throw new Exception("DO NOT CALL Triangle.fold")
  def unfold(index: Long): Vector[Long] = throw new Exception("DO NOT CALL Triangle.unfold")
  def indexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long], hints: Option[RangeComputeHints] = None): Seq[IndexRange] =
    throw new Exception("DO NOT CALL Triangle.indexRanges")

  def stackIndex(nextIndex: Long): Long = (index << 3) | (nextIndex & 7)

  def nextOrientation(transition: Int): Int = OrientationTransitions((orientation, transition))

  // returns the triangle (if available) corresponding to the next valid index
  def next: Option[Triangle] = {
    //println(s"Triangle.next from $bitString...")
    val transitions: Vector[Int] = bitString.sliding(3,3).map(bbb => BitStrings(bbb)).toVector
    var i = transitions.length - 1
    var result: Option[Triangle] = None
    while (i >= 0) {
      val t = transitions(i)
      if (i == 0) {
        if (t >= 7) return None
        val nextTransitions: Seq[Int] = Vector(t + 1) ++ Vector.fill(transitions.length - i - 1)(TransCenter)
        val nextIndex: Long = nextTransitions.foldLeft(0L)((acc, t) => {
          (acc << 3L) | t
        })
        return Option(invIndex(nextIndex, depth))
      } else {
        if (NextTransition.contains(t)) {
          val nextTransitions: Seq[Int] = (transitions.take(i) :+ NextTransition(t)) ++ Vector.fill(transitions.length - i - 1)(TransCenter)
          val nextIndex: Long = nextTransitions.foldLeft(0L)((acc, t) => {
            (acc << 3L) | t
          })
          return Option(invIndex(nextIndex, depth))
        }
      }
      i = i - 1
    }
    result
  }

  // the geo coordinates work funny for the region that includes the poles;
  // for a triangle tiling, this really just means the 8 portions (one per
  // oct face) that are furthest from the equator... in other words, no
  // matter your depth, these are the children you reached by only ever
  // recursing into the Apex starting from an oct face
  def isAllApex: Boolean = {
    if (depth == 1) return true
    for (i <- depth - 2 to 0 by -1) {
      val transition: Long = (index >> (3 * i)) & 3
      if (transition != TransApex) return false
    }
    true
  }

  /**
    * This is one of the key methods in the entire class; it defines how the geographic
    * and octahedral coordinate boundaries change as you recurse into child triangles.
    *
    * @param transition which child you are recursing into
    * @return the properly formatted child triangle
    */
  def child(transition: Int): Triangle = {
    val nextIndex = stackIndex(transition match {
      case TransCenter => Center
      case TransApex   => Corners(orientation)(0)
      case TransLL     => Corners(orientation)(2)
      case TransLR     => Corners(orientation)(1)
      case _ => throw new Exception("Invalid transition $transition")
    })

    val nextIsApexUp = isApexUp(nextOrientation(transition))

    //TODO figure out how to handle apex/base switching at the poles!
    val childBounds = transition match {
      case TransCenter =>
        //Extent(octX14, octX34), makeExtent(y0, yMid)
        TriBounds(bounds.octXBase.q2, bounds.octXBase.subH1, TriN.makeExtent(y0, yMid), nextIsApexUp)
      case TransApex =>
        //Extent(octX14, octX34), makeExtent(yMid, y1)
        TriBounds(bounds.octXBase.q2, bounds.octXBase.subH1, TriN.makeExtent(yMid, y1), nextIsApexUp)
      case TransLL =>
        //Extent(octX0, octXMid), makeExtent(y0, yMid)
        TriBounds(bounds.octXBase.q1, bounds.octXBase.subH0, TriN.makeExtent(y0, yMid), nextIsApexUp)
      case TransLR =>
        //Extent(octXMid, octX1), makeExtent(y0, yMid)
        TriBounds(bounds.octXBase.q3, bounds.octXBase.subH2, TriN.makeExtent(y0, yMid), nextIsApexUp)
      case _ =>
        throw new Exception(s"Invalid transition $transition")
    }

    transition match {
      case TransCenter => Triangle(nextIndex, nextOrientation(TransCenter), childBounds, depth + 1)
      case TransApex => Triangle(nextIndex, nextOrientation(TransApex), childBounds, depth + 1)
      case TransLL => Triangle(nextIndex, nextOrientation(TransLL), childBounds, depth + 1)
      case TransLR => Triangle(nextIndex, nextOrientation(TransLR), childBounds, depth + 1)
      case _ => throw new Exception(s"Invalid transition $transition")
    }
  }

  // assumes that you've ALREADY transformed X into octahedral form!
  def getTriangle(octX: Double, y: Double, maxDepth: Int): Triangle = {
    require(bounds.containsOctX(octX), s"X is out of bounds:  $octX notIn $bounds")
    require(bounds.containsY(y), s"Y is out of bounds:  $y notIn $bounds")

    if (maxDepth < 1) return this

    // apex
    if (Math.abs(y - y1) < Math.abs(y - y0)) {
      require(octX >= octX14 && octX <= octX34, s"Cannot recurse into apex, X $octX is out of range ($octX14, $octX34)")
      return child(TransApex).getTriangle(octX, y, maxDepth - 1)
    }

    // you know you're recursing into the base

    // check left and right
    if (octX < octXMid) {
      // check left
      if (octX < octX14) {
        return child(TransLL).getTriangle(octX, y, maxDepth - 1)
      }
      val xProbe: Double = mLeftInv * (y - y0) + octXMid
      if (xProbe >= octX0 && xProbe < octXMid && octX < xProbe) {
        return child(TransLL).getTriangle(octX, y, maxDepth - 1)
      }
    } else {
      // check right
      if (octX > octX34) {
        return child(TransLR).getTriangle(octX, y, maxDepth - 1)
      }
      val xProbe: Double = mRightInv * (y - y0) + octXMid
      if (xProbe >= octXMid && xProbe <= octX1 && octX >= xProbe) {
        return child(TransLR).getTriangle(octX, y, maxDepth - 1)
      }
    }

    // if you get this far, recurse center
    child(TransCenter).getTriangle(octX, y, maxDepth - 1)
  }

  def geoWkt: String = bounds.geoWkt(orientation)

  def octWkt: String = bounds.octWkt(orientation)
}

object TriN {
  val UseCompactIndex: Boolean = true

  val ApexBias: Double = 1e-14

  def min(x: Double*): Double = x.min
  def max(x: Double*): Double = x.max

  // PMLP for working with Extent[Double] ranges for geo-, oct-coordinates
  implicit class RichExtent(extent: Extent[Double]) {
    def span: Double = extent.max - extent.min

    def q0: Double = extent.min
    def q1: Double = 0.75 * extent.min + 0.25 * extent.max
    def q2: Double = 0.5 * (extent.min + extent.max)
    def q3: Double = 0.25 * extent.min + 0.75 * extent.max
    def q4: Double = extent.max

    val Alpha = 1e-10
    def insideMin: Double = q0 * Alpha + q1 * (1.0 - Alpha)
    def insideMax: Double = q4 * Alpha + q3 * (1.0 - Alpha)

    def subQ1: Extent[Double] = Extent(q0, q1)
    def subQ2: Extent[Double] = Extent(q1, q2)
    def subQ3: Extent[Double] = Extent(q2, q3)
    def subQ4: Extent[Double] = Extent(q3, q4)

    def subH0: Extent[Double] = Extent(q0, q2)
    def subH1: Extent[Double] = Extent(q1, q3)
    def subH2: Extent[Double] = Extent(q2, q4)

    def peq0: Extent[Double] = Extent(q0, q0)
    def peq1: Extent[Double] = Extent(q1, q1)
    def peq2: Extent[Double] = Extent(q2, q2)
    def peq3: Extent[Double] = Extent(q3, q3)
    def peq4: Extent[Double] = Extent(q4, q4)

    def apply(functor: Double => Double): Extent[Double] = Extent(functor(extent.min), functor(extent.max))
  }

  val TriTopLongitude: Longitude = Longitude(4)
  val TriTopLatitude: Latitude = Latitude(2)

  // corners
  val A = 1
  val B = 2
  val C = 4

  val Center = 0

  // possible triangle orientations...
  // apex up
  /*
  val A_BC = 0
  val B_CA = 1
  val C_AB = 2
  val A_CB = 3
  val B_AC = 4
  val C_BA = 5
  val AB_C = 6
  val BC_A = 7
  val CA_B = 8
  val AC_B = 9
  val BA_C = 10
  val CB_A = 11
  */

  // left hexagon
  val CB_A = 0
  val B_AC = 1
  val AC_B = 2
  val A_CB = 3
  val BA_C = 4
  val C_BA = 5
  // right hexagon
  val BC_A = 6
  val C_AB = 7
  val AB_C = 8
  val A_BC = 9
  val CA_B = 10
  val B_CA = 11

  val Orientations: Set[Int] = Set(
    A_BC, B_CA, C_AB, A_CB, B_AC, C_BA,
    AB_C, BC_A, CA_B, AC_B, BA_C, CB_A
  )

  val TransCenter = 0
  val TransApex = 1
  val TransLR = 2
  val TransLL = 4

  val Transitions: Set[Int] = Set(TransCenter, TransApex, TransLR, TransLL)

  val BitStrings: Map[String, Int] = (0 to 7).map(transition => {
    (transition.toBinaryString.reverse.padTo(3, "0").reverse.mkString(""), transition)
  }).toMap

  val NextTransition: Map[Int, Int] = Map(
    TransCenter -> TransApex,
    TransApex -> TransLR,
    TransLR -> TransLL
  )

  val MPos: Double = 2.0
  val MPosInv: Double = 1.0 / MPos
  val MNeg: Double = -2.0
  val MNegInv: Double = 1.0 / MNeg

  def makeExtent(a: Double, b: Double): Extent[Double] =
    if (a <= b) Extent(a, b)
    else Extent(b, a)

  def isValidOrientation(orientation: Int): Boolean = Orientations.contains(orientation)

  def isApexUp(orientation: Int): Boolean = orientation match {
    case A_BC | B_CA | C_AB | A_CB | B_AC | C_BA => true
    case AB_C | BC_A | CA_B | AC_B | BA_C | CB_A => false
    case _ => throw new Exception(s"Invalid orientation ($orientation)")
  }

  def H(orientation: Int): Int = orientation match {
    case A_BC => A_CB
    case B_CA => B_AC
    case C_AB => C_BA
    case A_CB => A_BC
    case B_AC => B_CA
    case C_BA => C_AB
    case AB_C => BA_C
    case BC_A => CB_A
    case CA_B => AC_B
    case AC_B => CA_B
    case BA_C => AB_C
    case CB_A => BC_A
    case _ => throw new Exception(s"Invalid orientation ($orientation)")
  }

  def V(orientation: Int): Int = orientation match {
    case A_BC => BC_A
    case B_CA => CA_B
    case C_AB => AB_C
    case A_CB => CB_A
    case B_AC => AC_B
    case C_BA => BA_C
    case AB_C => C_AB
    case BC_A => A_BC
    case CA_B => B_CA
    case AC_B => B_AC
    case BA_C => C_BA
    case CB_A => A_CB
    case _ => throw new Exception(s"Invalid orientation ($orientation)")
  }

  def CCW(orientation: Int): Int = orientation match {
    case A_BC => AB_C
    case B_CA => BC_A
    case C_AB => CA_B
    case A_CB => AC_B
    case B_AC => BA_C
    case C_BA => CB_A
    case AB_C => B_CA
    case BC_A => C_AB
    case CA_B => A_BC
    case AC_B => C_BA
    case BA_C => A_CB
    case CB_A => B_AC
    case _ => throw new Exception(s"Invalid orientation ($orientation)")
  }

  def CW(orientation: Int): Int = orientation match {
    case A_BC => CA_B
    case B_CA => AB_C
    case C_AB => BC_A
    case A_CB => BA_C
    case B_AC => CB_A
    case C_BA => AC_B
    case AB_C => A_BC
    case BC_A => B_CA
    case CA_B => C_AB
    case AC_B => A_CB
    case BA_C => B_AC
    case CB_A => C_BA
    case _ => throw new Exception(s"Invalid orientation ($orientation)")
  }

  val Corners: Map[Int, Vector[Int]] = Map(
    A_BC -> Vector(A, B, C),
    A_CB -> Vector(A, C, B),
    B_AC -> Vector(B, A, C),
    B_CA -> Vector(B, C, A),
    C_AB -> Vector(C, A, B),
    C_BA -> Vector(C, B, A),
    BC_A -> Vector(A, B, C),
    CB_A -> Vector(A, C, B),
    AC_B -> Vector(B, A, C),
    CA_B -> Vector(B, C, A),
    AB_C -> Vector(C, A, B),
    BA_C -> Vector(C, B, A)
  )

  val TransitionsByOrientationCorner: Map[(Int, Int), Int] = {
    val orientationCornersItr = CartesianProductIterable(Seq(Orientations.toSeq, Seq(Center, A, B, C)))
    orientationCornersItr.map {
      case Seq(orientation: Int, corner: Int) =>
        Corners(orientation.asInstanceOf[Int]) match {
          case Vector(c, _, _) if c == corner => (orientation.asInstanceOf[Int], corner.asInstanceOf[Int]) -> TransApex
          case Vector(_, c, _) if c == corner => (orientation.asInstanceOf[Int], corner.asInstanceOf[Int]) -> TransLR
          case Vector(_, _, c) if c == corner => (orientation.asInstanceOf[Int], corner.asInstanceOf[Int]) -> TransLL
          case _                              => (orientation.asInstanceOf[Int], corner.asInstanceOf[Int]) -> TransCenter
        }
    }
  }.toMap

  val OrientationTransitions: Map[(Int, Int), Int] = {
    val orientationTransitionsItr = CartesianProductIterable(Seq(Orientations.toSeq, Transitions.toSeq))
    orientationTransitionsItr.map {
      case Seq(orientation, transition) =>
        transition match {
          case TransApex => (orientation.asInstanceOf[Int], transition.asInstanceOf[Int]) -> H(orientation.asInstanceOf[Int])
          case TransCenter => (orientation.asInstanceOf[Int], transition.asInstanceOf[Int]) -> V(H(orientation.asInstanceOf[Int]))
          case TransLR => (orientation.asInstanceOf[Int], transition.asInstanceOf[Int]) -> CCW(CCW(H(orientation.asInstanceOf[Int])))
          case TransLL => (orientation.asInstanceOf[Int], transition.asInstanceOf[Int]) -> CW(CW(H(orientation.asInstanceOf[Int])))
          case _ => throw new Exception(s"Invalid transition in pair ($orientation, $transition)")
        }
    }.toMap
  }

  val OrientationNames: Map[Int, String] = Map(
    A_BC -> "A|BC",
    A_CB -> "A|CB",
    B_AC -> "B|AC",
    B_CA -> "B|CA",
    C_AB -> "C|AB",
    C_BA -> "C|BA",
    AB_C -> "AB|C",
    BA_C -> "BA|C",
    AC_B -> "AC|B",
    CA_B -> "CA|B",
    BC_A -> "BC|A",
    CB_A -> "CB|A"
  )

  def compactIndexLength(depth: Int): Int = Math.max(3, 3 + 2 * (depth - 1))

  def expandedIndex(index: Long, depth: Int): Long = expandedIndex(longToBinaryString(index, compactIndexLength(depth)))

  def expandedIndex(bits: String): Long = {
    val octFace: Long = java.lang.Long.parseLong(bits.take(3), 2)
    bits.drop(3).sliding(2, 2).foldLeft(octFace)((acc, bb) => {
      val bval: Long = java.lang.Long.parseLong(bb, 2)
      (acc << 3L) | (if (bval == 3) 4 else bval)
    })
  }

  def compactIndex(bitString: String): Long = {
    val longs: List[Long] = bitString.sliding(3, 3).map(bits => java.lang.Long.parseLong(bits, 2)).toList
    longs.tail.foldLeft(longs.head)((acc, bb) => (acc << 2) | (if (bb == 4) 3 else bb))
  }

  def topOrientation(ix: Int, iy: Int): Int = {
    require(ix >= 0 && ix < 4)
    require(iy >= 0 && iy < 2)
    (iy & 1, ix & 1) match {
      case (0, 0) => BC_A
      case (0, 1) => CB_A
      case (1, 0) => A_BC
      case (1, 1) => A_CB
    }
  }

  def createLowestIndex(depth: Int): Triangle = {
    invIndex(List.fill(depth)(TriN.TransCenter).foldLeft(0L)((acc, t) => (acc << 3L) | t), depth)
  }

  val OctahedronFaces: Vector[Triangle] = (0 to 7).map(face => invIndex(face, 1)).toVector

  // assume these are GEOGRAPHIC X-coordinates, /not/ octahedral X-coordinates
  def getRanges(geoXMin: Double, ymin: Double, geoXMax: Double, ymax: Double, maxDepth: Int): Seq[IndexRange] = {
    // translate the geo coordinates to octahedral coordinates
    val octXMin = Math.min(geoToOctX(ymin)(geoXMin), geoToOctX(ymax)(geoXMin))
    val octXMax = Math.max(geoToOctX(ymin)(geoXMax), geoToOctX(ymax)(geoXMax))

    // TODO remove after debugging!
    //println(f"TriN.getRanges($geoXMin%1.3f, $ymin%1.3f, $geoXMax%1.3f, $ymax%1.3f, $maxDepth%d)")
    //println(f"  ... becomes oct X $octXMin%1.3f, $octXMax%1.3f")

    // find the top-level octahedron faces that intersect the query region at all
    val faces: Vector[Triangle] = OctahedronFaces.filter(t => t.bounds.octOverlaps(octXMin, ymin, octXMax, ymax))
    if (faces.isEmpty) throw new Exception("Triangle query did not match any octahedral faces!")

    // recurse, and combine results
    faces.tail.foldLeft(faces.head.getRangesRecursively(octXMin, ymin, octXMax, ymax, maxDepth))((acc, face) => {
      val faceRanges: Seq[IndexRange] = face.getRangesRecursively(octXMin, ymin, octXMax, ymax, maxDepth)
      face.consolidateRanges((acc ++ faceRanges).iterator).toSeq
    })
  }

  // map the given (x, y) coordinate to one face of the octahedron
  // that wraps the entire (spherical) Earth
  // assume that this is a GEOGRAPHIC X-coordinate
  def getInitialTriangle(geoX: Double, y: Double): Triangle = {
    val ix = TriTopLongitude.index(Seq(geoX))
    val iy = TriTopLatitude.index(Seq(y))
    val orientation: Int = topOrientation(ix.toInt, iy.toInt)
    val idx: Long = (iy << 2) | (if (ix < 2) ix else 5 - ix)
    val x0: Double = -180.0 + 90.0 * ix
    val y0: Double = -90.0 + 90.0 * iy
    // the octahedral face is a weird case where the APEX bounds are also 90 degrees wide!
    val bounds: TriBounds = TriBounds(x0 + 45.0, Extent(x0, x0 + 90.0), Extent(y0, y0 + 90.0), isApexUp(orientation))
    Triangle(idx, orientation, bounds, 1)
  }

  def geoToOctX(y: Double)(geoX: Double): Double = {
    // the octahedral faces are all 90 degrees wide, so you just need to pick the
    // right 90-degree range that includes the given [geoX, geoX)
    val geoX0: Double = 90.0 * Math.floor((180.0 + Math.max(-180.0, Math.min(179.999999, geoX))) / 90.0) - 180.0

    // do stuff
    val alpha = Math.abs(y) / 90.0
    val xMid = geoX0 + 45.0
    val octX = alpha * xMid + (1.0 - alpha) * geoX

//    // given the range minimum, transform the geo-X into an oct-X
//    val pY = (90.0 - Math.abs(y)) / 90.0
//    val xMid = geoX0 + 45.0
//    val octX = xMid - pY * (xMid - geoX)

    octX
  }

  def octToGeoX(y: Double)(octX: Double): Double = {
    // the octahedral faces are all 90 degrees wide, so you just need to pick the
    // right 90-degree range that includes the given geoX
    val geoX0: Double = 90.0 * Math.floor((180.0 + Math.max(-180.0, Math.min(179.999999, octX))) / 90.0) - 180.0

    // do stuff
    val alpha = Math.abs(y) / 90.0
    val xMid = geoX0 + 45.0
    val geoX = if (alpha < 1.0) {
      (octX - alpha * xMid) / (1.0 - alpha)
    } else {
      // at the discontinuity, you have to go all the way to the extreme, so pick the
      // side that you favor (if any)
      if (octX <= xMid) geoX0 else geoX0 + 90.0
    }

//    // given the range minimum, transform the geo-X into an oct-X
//    val pY = (90.0 - Math.abs(y)) / 90.0
//    // mind the potential division-by-zero
//    if (pY < 1e-6) return octX
//    val xMid = geoX0 + 45.0
//    val geoX = xMid - (xMid - octX) / pY

    geoX
  }

  // assume you've given a GEOGRAPHIC X-coordinate
  def getTriangle(geoX: Double, y: Double, maxDepth: Int): Triangle = {
    require(maxDepth > 0, s"maxDepth must be at least one; found $maxDepth")
    val octFace = getInitialTriangle(geoX, y)
    if (maxDepth > 1) {
      val result = octFace.getTriangle(geoToOctX(y)(geoX), y, maxDepth - 1)
      result
    } else {
      val result = octFace
      result
    }
  }

  def longToBinaryString(long: Long, length: Int): String =
    long.toBinaryString.reverse.padTo(length, "0").reverse.mkString("")

  def indexBinaryString(index: Long, depth: Int): String = longToBinaryString(index, 3 * depth)

  def indexOctalString(index: Long, depth: Int): String = {
    (0 until depth).map(pos => ((index >> (3 * pos)) & 7L).toString).reverse.mkString("")
  }

  // assume you've given a GEOGRAPHIC X-coordinate
  def index(geoX: Double, y: Double, maxDepth: Int): Long = {
    getTriangle(geoX, y, maxDepth).index
  }

  def invIndex(idx: Long, depth: Int): Triangle = {
    require(depth > 0, s"Depth ($depth) must be at least 1")

    val octIdx = (idx >> (3 * (depth - 1))) & 7
    val octIdxX = octIdx & 3
    val octXi = if (octIdxX < 2) octIdxX else 5 - octIdxX
    val octYi = (octIdx >> 2) & 1
    val x0 = octXi * 90.0 - 180.0
    val y0 = octYi * 90.0 - 90.0
    val faceOrientation = topOrientation(octXi.toInt, octYi.toInt)
    val bounds: TriBounds = TriBounds(x0 + 45.0, Extent(x0, x0 + 90.0), Extent(y0, y0 + 90.0), isApexUp(faceOrientation))
    val octFace: Triangle = Triangle(octIdx, faceOrientation, bounds, 1)

    (depth - 2  to 0 by -1).foldLeft(octFace)((acc, i) => {
      val corner = (idx >> (3 * i)) & 7
      acc.child(TransitionsByOrientationCorner((acc.orientation, corner.toInt)))
    })
  }

  def invIndex(indexAsBinaryString: String): Triangle = {
    require(indexAsBinaryString != null)
    require(indexAsBinaryString.length % 3 == 0)
    val depth = indexAsBinaryString.length / 3
    val index = java.lang.Long.parseLong(indexAsBinaryString, 2)
    invIndex(index, depth)
  }

  // uses the "next" function within the Triangle class itself
  def iterator(depth: Int): Iterator[Triangle] = new Iterator[Triangle] {
    require(depth > 0, s"Depth ($depth) must be a positive integer")
    var triangle: Option[Triangle] = Option(createLowestIndex(depth))
    def hasNext: Boolean = triangle.isDefined
    def next(): Triangle = {
      val result = triangle.orNull
      triangle = triangle.get.next
      result
    }
  }

  def geoLineWkt(t0: Triangle, t1: Triangle): String =
    "LINESTRING(" +
      t0.bounds.geoXMid + " " + t0.bounds.yMid + ", " +
      t1.bounds.geoXMid + " " + t1.bounds.yMid +
      ")"

  def octLineWkt(t0: Triangle, t1: Triangle): String =
    "LINESTRING(" +
      t0.bounds.octXMid + " " + t0.bounds.yMid + ", " +
      t1.bounds.octXMid + " " + t1.bounds.yMid +
      ")"

  case class Degrees(degrees: Double) {
    def radians: Double = Math.toRadians(degrees)
  }

  case class Point(x: Degrees, y: Degrees) {
    def wkt: String = "POINT(" + x.degrees + " " + y.degrees + ")"
  }

}

object NamedLocations {
  case class LatLon(latitude: Double, longitude: Double)
  val LocationsByName: Map[String, LatLon] = Map(
    "Charlottesville" -> LatLon(38.0293, -78.4767),
    "CCRi" -> LatLon(38.075776, -78.495150),
    "Eichelberger" -> LatLon(38.054444, -78.688256),
    "Kunkel" -> LatLon(38.053707, -78.688751),
    "Price" -> LatLon(38.053373, -78.684150),
    "Uluru" -> LatLon(-25.344407, 131.036881),
    "Ipanema" -> LatLon(-22.986877, -43.208614),
    "Bogata" -> LatLon(4.592899, -74.123853),
    "Ravensburg" -> LatLon(47.776487, 9.610370),
    "McMurdo" -> LatLon(-77.840380, 166.685197),
    "Kilauea" -> LatLon(19.407757, -155.283448),
    "Reunion" -> LatLon(-21.097563, 55.480176),
    "Fukuoka" -> LatLon(33.593343, 130.351508),
    "Beijing" -> LatLon(39.908780, 116.397482),
    "Cabo San Lucas" -> LatLon(22.889986, -109.918799),
    "Hanga Roa" -> LatLon(-27.139735, -109.427335),
    "Wellington" -> LatLon(-41.288435, 174.777162),
    "South Pole" -> LatLon(-90.0, -179.0),
    "North Pole" -> LatLon(90.0 - 1e-10, 179.0),
    "Fiji" -> LatLon(-17.825785, 177.960057),
    "Tahiti" -> LatLon(-17.682831, -149.445502),
    "London" -> LatLon(51.501044, -0.124639),
    "Johannesburg" -> LatLon(-26.188086, 28.042677),
    "Riyadh" -> LatLon(24.698679, 46.672781),
    "Moscow" -> LatLon(55.752628, 37.623058),
    "Seattle" -> LatLon(47.622318, -122.316094),
    "Null Island" -> LatLon(0.0, 0.0)
  )
}

object TriTest extends App {
  import NamedLocations._

  var pw: PrintWriter = null

  // H must be invertible
  assert(A_BC == H(H(A_BC)), s"Uninvertible:  A_BC <> H(H(A_BC))")
  assert(B_CA == H(H(B_CA)), s"Uninvertible:  B_CA <> H(H(B_CA))")
  assert(C_AB == H(H(C_AB)), s"Uninvertible:  C_AB <> H(H(C_AB))")
  assert(A_CB == H(H(A_CB)), s"Uninvertible:  A_CB <> H(H(A_CB))")
  assert(B_AC == H(H(B_AC)), s"Uninvertible:  B_AC <> H(H(B_AC))")
  assert(C_BA == H(H(C_BA)), s"Uninvertible:  C_BA <> H(H(C_BA))")
  assert(AB_C == H(H(AB_C)), s"Uninvertible:  AB_C <> H(H(AB_C))")
  assert(BC_A == H(H(BC_A)), s"Uninvertible:  BC_A <> H(H(BC_A))")
  assert(CA_B == H(H(CA_B)), s"Uninvertible:  CA_B <> H(H(CA_B))")
  assert(AC_B == H(H(AC_B)), s"Uninvertible:  AC_B <> H(H(AC_B))")
  assert(BA_C == H(H(BA_C)), s"Uninvertible:  BA_C <> H(H(BA_C))")
  assert(CB_A == H(H(CB_A)), s"Uninvertible:  CB_A <> H(H(CB_A))")

  // V must be invertible
  assert(A_BC == V(V(A_BC)), s"Uninvertible:  A_BC <> V(V(A_BC))")
  assert(B_CA == V(V(B_CA)), s"Uninvertible:  B_CA <> V(V(B_CA))")
  assert(C_AB == V(V(C_AB)), s"Uninvertible:  C_AB <> V(V(C_AB))")
  assert(A_CB == V(V(A_CB)), s"Uninvertible:  A_CB <> V(V(A_CB))")
  assert(B_AC == V(V(B_AC)), s"Uninvertible:  B_AC <> V(V(B_AC))")
  assert(C_BA == V(V(C_BA)), s"Uninvertible:  C_BA <> V(V(C_BA))")
  assert(AB_C == V(V(AB_C)), s"Uninvertible:  AB_C <> V(V(AB_C))")
  assert(BC_A == V(V(BC_A)), s"Uninvertible:  BC_A <> V(V(BC_A))")
  assert(CA_B == V(V(CA_B)), s"Uninvertible:  CA_B <> V(V(CA_B))")
  assert(AC_B == V(V(AC_B)), s"Uninvertible:  AC_B <> V(V(AC_B))")
  assert(BA_C == V(V(BA_C)), s"Uninvertible:  BA_C <> V(V(BA_C))")
  assert(CB_A == V(V(CB_A)), s"Uninvertible:  CB_A <> V(V(CB_A))")

  // CCW must go all the way around
  assert(A_BC == CCW(CCW(CCW(CCW(CCW(CCW(A_BC)))))), s"Uninvertible:  A_BC <> CCW(CCW(CCW(CCW(CCW(CCW(A_BC))))))")
  assert(B_CA == CCW(CCW(CCW(CCW(CCW(CCW(B_CA)))))), s"Uninvertible:  B_CA <> CCW(CCW(CCW(CCW(CCW(CCW(B_CA))))))")
  assert(C_AB == CCW(CCW(CCW(CCW(CCW(CCW(C_AB)))))), s"Uninvertible:  C_AB <> CCW(CCW(CCW(CCW(CCW(CCW(C_AB))))))")
  assert(A_CB == CCW(CCW(CCW(CCW(CCW(CCW(A_CB)))))), s"Uninvertible:  A_CB <> CCW(CCW(CCW(CCW(CCW(CCW(A_CB))))))")
  assert(B_AC == CCW(CCW(CCW(CCW(CCW(CCW(B_AC)))))), s"Uninvertible:  B_AC <> CCW(CCW(CCW(CCW(CCW(CCW(B_AC))))))")
  assert(C_BA == CCW(CCW(CCW(CCW(CCW(CCW(C_BA)))))), s"Uninvertible:  C_BA <> CCW(CCW(CCW(CCW(CCW(CCW(C_BA))))))")
  assert(AB_C == CCW(CCW(CCW(CCW(CCW(CCW(AB_C)))))), s"Uninvertible:  AB_C <> CCW(CCW(CCW(CCW(CCW(CCW(AB_C))))))")
  assert(BC_A == CCW(CCW(CCW(CCW(CCW(CCW(BC_A)))))), s"Uninvertible:  BC_A <> CCW(CCW(CCW(CCW(CCW(CCW(BC_A))))))")
  assert(CA_B == CCW(CCW(CCW(CCW(CCW(CCW(CA_B)))))), s"Uninvertible:  CA_B <> CCW(CCW(CCW(CCW(CCW(CCW(CA_B))))))")
  assert(AC_B == CCW(CCW(CCW(CCW(CCW(CCW(AC_B)))))), s"Uninvertible:  AC_B <> CCW(CCW(CCW(CCW(CCW(CCW(AC_B))))))")
  assert(BA_C == CCW(CCW(CCW(CCW(CCW(CCW(BA_C)))))), s"Uninvertible:  BA_C <> CCW(CCW(CCW(CCW(CCW(CCW(BA_C))))))")
  assert(CB_A == CCW(CCW(CCW(CCW(CCW(CCW(CB_A)))))), s"Uninvertible:  CB_A <> CCW(CCW(CCW(CCW(CCW(CCW(CB_A))))))")

  // CW must go all the way around
  assert(A_BC == CW(CW(CW(CW(CW(CW(A_BC)))))), s"Uninvertible:  A_BC <> CW(CW(CW(CW(CW(CW(A_BC))))))")
  assert(B_CA == CW(CW(CW(CW(CW(CW(B_CA)))))), s"Uninvertible:  B_CA <> CW(CW(CW(CW(CW(CW(B_CA))))))")
  assert(C_AB == CW(CW(CW(CW(CW(CW(C_AB)))))), s"Uninvertible:  C_AB <> CW(CW(CW(CW(CW(CW(C_AB))))))")
  assert(A_CB == CW(CW(CW(CW(CW(CW(A_CB)))))), s"Uninvertible:  A_CB <> CW(CW(CW(CW(CW(CW(A_CB))))))")
  assert(B_AC == CW(CW(CW(CW(CW(CW(B_AC)))))), s"Uninvertible:  B_AC <> CW(CW(CW(CW(CW(CW(B_AC))))))")
  assert(C_BA == CW(CW(CW(CW(CW(CW(C_BA)))))), s"Uninvertible:  C_BA <> CW(CW(CW(CW(CW(CW(C_BA))))))")
  assert(AB_C == CW(CW(CW(CW(CW(CW(AB_C)))))), s"Uninvertible:  AB_C <> CW(CW(CW(CW(CW(CW(AB_C))))))")
  assert(BC_A == CW(CW(CW(CW(CW(CW(BC_A)))))), s"Uninvertible:  BC_A <> CW(CW(CW(CW(CW(CW(BC_A))))))")
  assert(CA_B == CW(CW(CW(CW(CW(CW(CA_B)))))), s"Uninvertible:  CA_B <> CW(CW(CW(CW(CW(CW(CA_B))))))")
  assert(AC_B == CW(CW(CW(CW(CW(CW(AC_B)))))), s"Uninvertible:  AC_B <> CW(CW(CW(CW(CW(CW(AC_B))))))")
  assert(BA_C == CW(CW(CW(CW(CW(CW(BA_C)))))), s"Uninvertible:  BA_C <> CW(CW(CW(CW(CW(CW(BA_C))))))")
  assert(CB_A == CW(CW(CW(CW(CW(CW(CB_A)))))), s"Uninvertible:  CB_A <> CW(CW(CW(CW(CW(CW(CB_A))))))")

  // CW must be the inverse of CCW
  assert(A_BC == CW(CCW(A_BC)), s"Uninvertible:  A_BC <> CW(CCW(A_BC))")
  assert(B_CA == CW(CCW(B_CA)), s"Uninvertible:  B_CA <> CW(CCW(B_CA))")
  assert(C_AB == CW(CCW(C_AB)), s"Uninvertible:  C_AB <> CW(CCW(C_AB))")
  assert(A_CB == CW(CCW(A_CB)), s"Uninvertible:  A_CB <> CW(CCW(A_CB))")
  assert(B_AC == CW(CCW(B_AC)), s"Uninvertible:  B_AC <> CW(CCW(B_AC))")
  assert(C_BA == CW(CCW(C_BA)), s"Uninvertible:  C_BA <> CW(CCW(C_BA))")
  assert(AB_C == CW(CCW(AB_C)), s"Uninvertible:  AB_C <> CW(CCW(AB_C))")
  assert(BC_A == CW(CCW(BC_A)), s"Uninvertible:  BC_A <> CW(CCW(BC_A))")
  assert(CA_B == CW(CCW(CA_B)), s"Uninvertible:  CA_B <> CW(CCW(CA_B))")
  assert(AC_B == CW(CCW(AC_B)), s"Uninvertible:  AC_B <> CW(CCW(AC_B))")
  assert(BA_C == CW(CCW(BA_C)), s"Uninvertible:  BA_C <> CW(CCW(BA_C))")
  assert(CB_A == CW(CCW(CB_A)), s"Uninvertible:  CB_A <> CW(CCW(CB_A))")

  // CCW must be the inverse of CW
  assert(A_BC == CCW(CW(A_BC)), s"Uninvertible:  A_BC <> CCW(CW(A_BC))")
  assert(B_CA == CCW(CW(B_CA)), s"Uninvertible:  B_CA <> CCW(CW(B_CA))")
  assert(C_AB == CCW(CW(C_AB)), s"Uninvertible:  C_AB <> CCW(CW(C_AB))")
  assert(A_CB == CCW(CW(A_CB)), s"Uninvertible:  A_CB <> CCW(CW(A_CB))")
  assert(B_AC == CCW(CW(B_AC)), s"Uninvertible:  B_AC <> CCW(CW(B_AC))")
  assert(C_BA == CCW(CW(C_BA)), s"Uninvertible:  C_BA <> CCW(CW(C_BA))")
  assert(AB_C == CCW(CW(AB_C)), s"Uninvertible:  AB_C <> CCW(CW(AB_C))")
  assert(BC_A == CCW(CW(BC_A)), s"Uninvertible:  BC_A <> CCW(CW(BC_A))")
  assert(CA_B == CCW(CW(CA_B)), s"Uninvertible:  CA_B <> CCW(CW(CA_B))")
  assert(AC_B == CCW(CW(AC_B)), s"Uninvertible:  AC_B <> CCW(CW(AC_B))")
  assert(BA_C == CCW(CW(BA_C)), s"Uninvertible:  BA_C <> CCW(CW(BA_C))")
  assert(CB_A == CCW(CW(CB_A)), s"Uninvertible:  CB_A <> CCW(CW(CB_A))")

  // test the geo-to-oct conversion both ways
  for (geoX <- -90.0 to -1e-4 by Math.PI; y <- 0.0 to 90.0 - 1e-4 by Math.PI) {
    val octX = geoToOctX(y)(geoX)
    val geoX2 = octToGeoX(y)(octX)
    //println(s"Y $y:  geoX $geoX -> octX $octX -> geoX $geoX2:  difference is ${Math.abs(geoX - geoX2)}")
    assert(Math.abs(geoX - geoX2) <= 1e-10, s"Y $y for geoX $geoX -> octX $octX -> geoX $geoX2:  difference is ${Math.abs(geoX - geoX2)}")
  }

  // index and inverse-index must be -- wait for it... -- inverses
  for (depth <- 1 to 21; loc <- LocationsByName) {
    loc match {
      case (name, LatLon(y, x)) =>
        val tFwd = getTriangle(x, y, depth)
        val idx = index(x, y, depth)
        assert(tFwd.index == idx, s"Mismatched indexes; triangle ${tFwd.index}, direct $idx")
        val tRev = invIndex(idx, depth)
        assert(tRev.index == idx, s"Mismatched reverse indexes; triangle ${tRev.index}, direct $idx")
    }
  }

  // ensure that you can encode many different triangles, and that the geo/oct conversion works round-trip
  def testInitialTri(n: Int, geoX: Double, y: Double): Unit = {
    val octFace = getTriangle(geoX, y, 1)
    val bIndex = octFace.index.toBinaryString.reverse.padTo(3, "0").reverse.mkString("")
    val octX = geoToOctX(y)(geoX)
    val x2 = octToGeoX(y)(octX)
    //pw.println(f"$n%d\tPOINT($geoX%1.4f $y%1.4f)\tPOINT($octX%1.4f $y%1.4f)\t$bIndex%s")
    if (Math.abs(y) < 89.5) {
      require(Math.abs(x2 - geoX) <= 1e-6, f"Failed to satisfy XT inverse:  $geoX%1.6f <> $octX%1.6f")
    }
  }
  val xs = (-180.0 to 179.999 by 22.5).zipWithIndex
  val ys = (-90.0 to 89.999 by 11.25*0.5).zipWithIndex
  val ny = ys.length
  for (xx <- xs; yy <- ys) {
    testInitialTri(yy._2 * ny + xx._2, xx._1, yy._1)
  }

  // drive "isInitialApex"
  for (oct <- 0 to 7) assert(invIndex(indexBinaryString(oct, 1)).isAllApex)
  for (oct <- 0 to 7; trans0 <- Transitions) {
    val t = invIndex(indexBinaryString(oct << 3 | trans0, 2)).isAllApex
    assert((trans0 != TransApex) ^ t)
  }
  for (oct <- 0 to 7; trans0 <- Transitions; trans1 <- Transitions) {
    val t = invIndex(indexBinaryString((oct << 3 | trans0) << 3 | trans1, 3)).isAllApex
    assert((trans0 != TransApex || trans1 != TransApex) ^ t)
  }
  for (oct <- 0 to 7; trans0 <- Transitions; trans1 <- Transitions; trans2 <- Transitions) {
    val t = invIndex(indexBinaryString(((oct << 3 | trans0) << 3 | trans1) << 3 | trans2, 4)).isAllApex
    assert((trans0 != TransApex || trans1 != TransApex || trans2 != TransApex) ^ t)
  }

  def testOctOverlap(t: Triangle, octX: Double, y: Double, expectation: Boolean): Unit = {
    val r = Rectangle(Extent(octX, octX, incMin = true, incMax = true), Extent(y, y, incMin = true, incMax = true))
    println(s"  test overlap expects $expectation\n    r:  $r\n    t:  $t\n    o:  ${t.bounds.octOverlaps(r)}")
    assert(t.bounds.octOverlaps(r) == expectation, "Triangle.bounds.octOverlaps(Rectangle) did not meet expectations.")
  }

  // test the "triangle overlaps (oct) rectangle" functions
  {
    val target: LatLon = LocationsByName("Eichelberger")
    val t: Triangle = TriN.getTriangle(target.longitude, target.latitude, 4)
    println(s"Test triangle:  $t")
    // obviously wrong
    testOctOverlap(t, -1.0, -1.0, expectation = false)
    testOctOverlap(t, -1.0, 1.0, expectation = false)
    testOctOverlap(t, 1.0, -1.0, expectation = false)
    testOctOverlap(t, 1.0, 1.0, expectation = false)
    // less obviously wrong
    testOctOverlap(t, -65.0, 44.9, expectation = false)
    testOctOverlap(t, -73.0, 44.0, expectation = false)
    // obviously right
    testOctOverlap(t, -67.5, 40.0, expectation = true)
    testOctOverlap(t, -68.0, 40.0, expectation = true)
    // less obviously right
    testOctOverlap(t, -62.0, 34.0, expectation = true)
    testOctOverlap(t, -73.0, 34.0, expectation = true)
  }

  // additional triangle-overlaps tests observed to fail during use
  {
    val depth: Int = 8
    val point = Point(Degrees(-78.688256), Degrees(38.054444))
    val geoX0: Double = Math.floor(point.x.degrees)
    val geoX1: Double = Math.ceil(point.x.degrees)
    val y0: Double = Math.floor(point.y.degrees)
    val y1: Double = Math.ceil(point.y.degrees)
    val octX0: Double = Math.min(geoToOctX(y0)(geoX0), geoToOctX(y1)(geoX0))
    val octX1: Double = Math.max(geoToOctX(y0)(geoX1), geoToOctX(y1)(geoX1))
    val r = Rectangle(Extent(octX0, octX1, incMin = true, incMax = true), Extent(y0, y1, incMin = true, incMax = true))
    val g2X0: Double = Math.min(TriN.octToGeoX(y0)(octX0), TriN.octToGeoX(y1)(octX0))
    val g2X1: Double = Math.max(TriN.octToGeoX(y0)(octX1), TriN.octToGeoX(y1)(octX1))
    ;
    val points: Seq[Point] = Seq(
      Point(Degrees(g2X1), Degrees(y0))
    )
    for (point <- points) {
      println(s"\nTest triangle...")
      println(s"  point")
      println(s"    geo:  X ${point.x.degrees}, Y ${point.y.degrees}")
      println(s"    oct:  X ${TriN.geoToOctX(point.y.degrees)(point.x.degrees)}, Y ${point.y.degrees}")
      val t: Triangle = TriN.getTriangle(point.x.degrees, point.y.degrees, depth)
      println(s"  triangle $t")
      println(s"  is apex up?  ${TriN.isApexUp(t.orientation)}")
      println(s"  rectangle $r")
      assert(t.bounds.octOverlaps(r), "Triangle should have intersected rectangle, but did not")
    }
  }

  // ensure that you can convert between compact and expanded index forms
  for (depth <- Seq(5); t <- iterator(depth)) {
    val nativeIndex = t.index
    val nativeIndexString = TriN.indexBinaryString(nativeIndex, depth)
    val nativeIndexLength = 3 * depth
    val compactIndex = TriN.compactIndex(nativeIndexString)
    val compactIndexLength = 3 + (2 * (depth - 1))
    val compactIndexString = TriN.longToBinaryString(compactIndex, compactIndexLength)
    val expandedIndex = TriN.expandedIndex(compactIndexString)
    val expandedIndexString = TriN.indexBinaryString(expandedIndex, depth)
    require(nativeIndexString.length == nativeIndexLength, s"Expected native index string to have $nativeIndexLength character(s), found ${nativeIndexString.length} instead")
    require(compactIndexString.length == compactIndexLength, s"Expected compact index string to have $compactIndexLength character(s), found ${compactIndexString.length} instead")
    require(expandedIndexString.length == nativeIndexLength, s"Expected expanded index string to have $nativeIndexLength character(s), found ${expandedIndexString.length} instead")
    //println(s"Native $nativeIndex $nativeIndexString -> $compactIndex $compactIndexString -> $expandedIndex $expandedIndexString")
    require(expandedIndex == nativeIndex, s"Native index ($nativeIndex, $nativeIndexString) != expanded index ($expandedIndex, $expandedIndexString)")
    val tInv = TriN.invIndex(expandedIndex, depth)
    require(t.index == tInv.index, s"Round-trip index failed:  $t != $tInv")
  }

  // test query known to return 8 triangles if working correctly
  {
    val depth = 8
    val point = Point(Degrees(-78.688256), Degrees(38.054444))
    val geoX0: Double = Math.floor(point.x.degrees)
    val geoX1: Double = Math.ceil(point.x.degrees)
    val y0: Double = Math.floor(point.y.degrees)
    val y1: Double = Math.ceil(point.y.degrees)
    val octX0: Double = Math.min(geoToOctX(y0)(geoX0), geoToOctX(y1)(geoX0))
    val octX1: Double = Math.max(geoToOctX(y0)(geoX1), geoToOctX(y1)(geoX1))
    val tOctFace = TriN.getTriangle(point.x.degrees, point.y.degrees, 1)
    val ranges = tOctFace.getRangesRecursively(octX0, y0, octX1, y1, depth)
    val numTriangles = ranges.map(range => range.upper - range.lower + 1).sum
    require(numTriangles == 8, s"Query for geo ($geoX0, $y0, $geoX1, $y1) expected 8 triangles, but got $numTriangles")
  }

  // another real case that failed
  /*
    Processing sample 0:  POINT(-20.865941382897603 85.04820360224821), POINT(-19.865941382897603 86.04820360224821)...

    Exception in thread "main" java.lang.Exception: Parent
      Triangle(5522195474,1,TriBounds(-43.59375,Extent(-43.6376953125,-43.5498046875,true,false),Extent(84.990234375,85.078125,true,false),true),11)
    overlaps oct bounds
      (-43.94030126770893, 85.04820360224821, -43.617125100877224, 86.04820360224821)
    but no child does:
      Triangle(44177563792,10,TriBounds(-43.59375,Extent(-43.61572265625,-43.57177734375,true,false),Extent(84.990234375,85.0341796875,true,false),false),12)
      Triangle(44177563794,11,TriBounds(-43.59375,Extent(-43.61572265625,-43.57177734375,true,false),Extent(85.0341796875,85.078125,true,false),true),12)
      Triangle(44177563793,7,TriBounds(-43.57177734375,Extent(-43.59375,-43.5498046875,true,false),Extent(84.990234375,85.0341796875,true,false),true),12)
      Triangle(44177563796,9,TriBounds(-43.61572265625,Extent(-43.6376953125,-43.59375,true,false),Extent(84.990234375,85.0341796875,true,false),true),12)
  {
    val t = TriN.invIndex(5522195474L, 11)
    val tApex = TriN.invIndex(44177563794L, 12)
    println("\n\n" + "="*80)
    println(s"\nYATT:  $t")
    val octXMin = -43.94030126770893
    val yMin = 85.04820360224821
    val octXMax = -43.617125100877224
    val yMax = 86.04820360224821
    println("computing parent overlap...")
    println(s"  ${t.bounds.octOverlaps(octXMin, yMin, octXMax, yMax)}")
    println("computing child overlap...")
    println(s"  ${tApex.bounds.octOverlaps(octXMin, yMin, octXMax, yMax)}")
    val ranges = t.getRangesRecursively(octXMin, yMin, octXMax, yMax, 12)
    println(s"ranges:  ${ranges.map(range => range.toString)}")
  }
  */

  def poly(geoMinX: Double, minY: Double, geoMaxX: Double, maxY: Double): String =
    "POLYGON((" +
    s"$geoMinX $minY, " +
    s"$geoMinX $maxY, " +
    s"$geoMaxX $maxY, " +
    s"$geoMaxX $minY, " +
    s"$geoMinX $minY" +
    "))"

  // dump files out that we want to visualize
  // (these are not unit tests so much as manual validation tests)
  try {
    // write out a query and the resulting ranges ask WKTs
    for (depth <- Seq(8, 9, 12)) {
      val point = Point(Degrees(-78.688256), Degrees(38.054444))
      val geoX0: Double = Math.floor(point.x.degrees)
      val geoX1: Double = Math.ceil(point.x.degrees)
      val y0: Double = Math.floor(point.y.degrees)
      val y1: Double = Math.ceil(point.y.degrees)
      val octX0: Double = Math.min(geoToOctX(y0)(geoX0), geoToOctX(y1)(geoX0))
      val octX1: Double = Math.max(geoToOctX(y0)(geoX1), geoToOctX(y1)(geoX1))
      val tOctFace = TriN.getTriangle(point.x.degrees, point.y.degrees, 1)
      println(s"depth $depth, point $point, oct face $tOctFace")
      pw = new PrintWriter(new FileWriter(s"test-query-d${depth}.txt"))
      pw.println(s"nature\tindex\tbits\tgeo_wkt\toct_wkt")
      pw.println(s"query\t\t\t${poly(geoX0, y0, geoX1, y1)}\t${poly(octX0, y0, octX1, y1)}")
      val ranges = tOctFace.getRangesRecursively(octX0, y0, octX1, y1, depth)
      for (range <- ranges; rangeIndex <- range.lower to range.upper) {
        val tIndex = TriN.expandedIndex(rangeIndex, depth)
        val t = TriN.invIndex(tIndex, depth)
        //println(s"range index $rangeIndex, expanded $tIndex, triangle $t")
        pw.println(s"range\t${t.index}\t${t.bitString}\t\t${t.octWkt}")
      }
      pw.close()
    }

    pw = new PrintWriter(new FileWriter("test-points.txt"))
    pw.println("name\tgeo_wkt\toct_wkt")
    LocationsByName.foreach {
      case (name, LatLon(y, geoX)) =>
        val t = getInitialTriangle(geoX, y)
        pw.println(s"$name\tPOINT($geoX $y)\tPOINT(${geoToOctX(y)(geoX)} $y)")
    }
    pw.close()

    pw = new PrintWriter(new FileWriter("test-triangles.txt"))
    pw.println("index\toct_wkt\tgeo_wkt")
    for (row <- 0 to 1; col <- 0 to 3) {
      val x0 = -180.0 + 90.0 * col
      val y0 = -90.0 + 90.0 * row
      val x1 = x0 + 90.0 - 1e-6
      val y1 = y0 + 90.0 - 1e-6
      val xMid = 0.5 * (x0 + x1)
      val yMid = 0.5 * (y0 + y1)
      val t = getInitialTriangle(xMid, yMid)
      pw.println(s"${t.bitString}\t${t.octWkt}\t${t.geoWkt}")
    }
    pw.close()

    pw = new PrintWriter(new FileWriter("test-index.txt"))
    pw.println("depth\tindex_dec\tindex_bits\toct_wkt\tgeo_wkt\tarea_deg2\tarea_m2")
    val target: LatLon = LocationsByName("Eichelberger")
    (1 to 21).foldLeft((target.longitude, target.latitude))((acc, depth) => acc match {
      case (x, y) =>
        val t = getTriangle(target.longitude, target.latitude, depth)
        pw.println(depth + "\t" + t.index + "\t" + t.bitString + "\t" + t.octWkt + "\t" + t.geoWkt + "\t" + t.bounds.areaSquareDegrees + "\t" + t.bounds.areaSquareMeters)
        (x, y)
    })
    pw.close()

    pw = new PrintWriter(new FileWriter("test-tiles.txt"))
    pw.println("orientation\tindex_dec\tindex_bits\toct_wkt\tgeo_wkt\tis_all_apex")
    for (t <- iterator(3)) {
      pw.println(s"${OrientationNames(t.orientation)}\t${t.index}\t${indexOctalString(t.index, t.depth)}\t${t.octWkt}\t${t.geoWkt}\t${if (t.isAllApex) 1 else 0}")
    }
    pw.close()

    pw = new PrintWriter(new FileWriter("test-progression.txt"))
    pw.println("index_from\tindex_to\toct_wkt\tgeo_wkt")
    for (t <- iterator(4).sliding(2, 1)) t match {
      case Seq(t0, t1) =>
        pw.println(s"${t0.index}\t${t1.index}\t${geoLineWkt(t0, t1)}\t${octLineWkt(t0, t1)}")
    }
    pw.close()

    // validate the 2D geo/oct coordinate translation
    pw = new PrintWriter(new FileWriter("test-oct2geo.txt"))
    pw.println("depth\tindex_oct\toct_wkt\tgeo_wkt")
    val target2: LatLon = LocationsByName("Eichelberger")
    (1 to 3).foldLeft((target2.longitude, target2.latitude))((acc, depth) => acc match {
      case (x, y) =>
        val t = getTriangle(target2.longitude, target2.latitude, depth)
        pw.println(depth + "\t" + indexOctalString(t.index, depth) + "\t" + t.octWkt + "\t" + t.geoWkt)
        (x, y)
    })
    pw.close()

  } finally {
    pw.close()
  }
}