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

// although this will often have a degenerate Apex Geo X (when the shape is a proper triangle),
// it can also describe a real rectangle (as for an octahedron face) as well as a trapezoid
case class TriBounds(geoXApex: Extent[Double], geoXBase: Extent[Double], octX: Extent[Double], Y: Extent[Double]) {
  val octX0: Double = octX.min
  val octX1: Double = octX.max
  val octXMid: Double = octX.q2
  val octX14: Double = octX.q1
  val octX34: Double = octX.q3

  val geoX: Extent[Double] = Extent(Math.min(geoXApex.min, geoXBase.min), Math.max(geoXApex.max, geoXBase.max))
  val geoX0: Double = geoX.min
  val geoX1: Double = geoX.max
  val geoXMid: Double = geoX.q2
  val geoX14: Double = geoX.q1
  val geoX34: Double = geoX.q3

  val y0: Double = Y.min
  val y1: Double = Y.max
  val yMid: Double = Y.q2

  def geoCenter: (Double, Double) = (geoXMid, yMid)

  def octCenter: (Double, Double) = (octXMid, yMid)

  def geoXAtMidY: Extent[Double] = Extent(
    0.5 * (geoXApex.min + geoXBase.min),
    0.5 * (geoXApex.max + geoXBase.max)
  )

  // this does allow for error; consider the case where a (proper, geo) rectangle
  // has a lower-right corner that is almost (but not) touching the left ascent of
  // this trapezoid:  this routine will falsely claim that the shapes everlap
  def geoOverlaps(xmin: Double, ymin: Double, xmax: Double, ymax: Double): Boolean = {
    // coarse definition
    if (geoX0 > xmax) return false
    if (geoX1 < xmin) return false
    if (y0 > ymax) return false
    if (y1 < ymin) return false
    true
  }

  def containsGeoPoint(x: Double, y: Double): Boolean =
    geoOverlaps(x, y, x, y)

  def containsGeoX(x: Double): Boolean = x >= geoX0 && x <= geoX1

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
      s"POLYGON((${octXMid} ${Y.max}, ${octX1} ${Y.min}, ${octX0} ${Y.min}, ${octXMid} ${Y.max}))"
    case AB_C | AC_B | BA_C | BC_A | CA_B | CB_A =>
      // apex on bottom
      s"POLYGON((${octX0} ${Y.max}, ${octX1} ${Y.max}, ${octXMid} ${Y.min}, ${octX0} ${Y.max}))"
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
  * @param index
  * @param orientation
  * @param bounds:  the geographic and octahedral bounds of this triangle
  * @param depth
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

  val name: String = "Triangle"

  def bitString: String = indexBinaryString(index, depth)

  // degenerate for the Curve contract
  override def toIndex(x: Double, y: Double): Long = {
    TriN.index(x, y, depth)
  }

  // degenerate for the Curve contract
  override def toPoint(index: Long): (Double, Double) = bounds.geoCenter

  // every triangle has four smaller sub-triangle children
  def childTriangles: Seq[Triangle] = Transitions.toSeq.map(child)

  // called by TriN and itself
  def getRangesRecursively(xmin: Double, ymin: Double, xmax: Double, ymax: Double, maxDepth: Int): Seq[IndexRange] = {
    // check stop conditions
    if (depth == maxDepth) {
      // if you've gotten here, then your (single) index is the only thing to return
      return Seq(IndexRange(index, index, contained = true))
    }
    if (depth > maxDepth) throw new Exception("Should not be able to recurse this far")

    // find your children that intersect the query area
    val subs = childTriangles.filter(_.bounds.geoOverlaps(xmin, ymin, xmax, ymax))
    if (subs.isEmpty) {
      throw new Exception(s"Parent overlaps geo bounds, but no child does")
      //return Seq[IndexRange]()
    }

    // recurse, and combine the results
    subs.tail.foldLeft(subs.head.getRangesRecursively(xmin, ymin, xmax, ymax, maxDepth))((acc, t) => {
      t.consolidateRanges((acc ++ t.getRangesRecursively(xmin, ymin, xmax, ymax, maxDepth)).iterator).toSeq
    })
  }

  // degenerate for the Curve contract
  override def toRanges(xmin: Double, ymin: Double, xmax: Double, ymax: Double, hints: Option[RangeComputeHints] = None): Seq[IndexRange] = {
    require(xmin <= xmax)
    require(ymin <= ymax)

    // defer to the companion object to do this planning top-down
    TriN.getRanges(xmin, ymin, xmax, ymax, depth)
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

    //TODO figure out how to handle apex/base switching at the poles!
    val childBounds = transition match {
      case TransCenter =>
        //Extent(octX14, octX34), makeExtent(y0, yMid)
        if (isAllApex) {
          TriBounds(bounds.geoXBase, bounds.geoXBase, bounds.octX.subH1, TriN.makeExtent(y0, yMid))
        } else {
          TriBounds(bounds.geoXBase.peq2, bounds.geoXAtMidY, bounds.octX.subH1, TriN.makeExtent(y0, yMid))
        }
      case TransApex =>
        //Extent(octX14, octX34), makeExtent(yMid, y1)
        if (isAllApex) {
          TriBounds(bounds.geoXApex, bounds.geoXBase, bounds.octX.subH1, TriN.makeExtent(yMid, y1))
        } else {
          TriBounds(bounds.geoXApex, bounds.geoXAtMidY, bounds.octX.subH1, TriN.makeExtent(yMid, y1))
        }
      case TransLL =>
        //Extent(octX0, octXMid), makeExtent(y0, yMid)
        if (isAllApex) {
          TriBounds(bounds.geoXAtMidY.peq0, bounds.geoXBase.subH0, bounds.octX.subH0, TriN.makeExtent(y0, yMid))
        } else {
          TriBounds(bounds.geoXAtMidY.peq0, bounds.geoXAtMidY.subH0, bounds.octX.subH0, TriN.makeExtent(y0, yMid))
        }
      case TransLR =>
        //Extent(octXMid, octX1), makeExtent(y0, yMid)
        if (isAllApex) {
          TriBounds(bounds.geoXAtMidY.peq4, bounds.geoXBase.subH2, bounds.octX.subH2, TriN.makeExtent(y0, yMid))
        } else {
          TriBounds(bounds.geoXAtMidY.peq4, bounds.geoXAtMidY.subH2, bounds.octX.subH2, TriN.makeExtent(y0, yMid))
        }
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
  // PMLP for working with Extent[Double] ranges for geo-, oct-coordinates
  implicit class RichExtent(extent: Extent[Double]) {
    def span: Double = extent.max - extent.min

    def q0: Double = extent.min
    def q1: Double = 0.75 * extent.min + 0.25 * extent.max
    def q2: Double = 0.5 * (extent.min + extent.max)
    def q3: Double = 0.25 * extent.min + 0.75 * extent.max
    def q4: Double = extent.max

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
  def getRanges(xmin: Double, ymin: Double, xmax: Double, ymax: Double, maxDepth: Int): Seq[IndexRange] = {
    // find the top-level octahedron faces that intersect the query region at all
    val faces: Vector[Triangle] = OctahedronFaces.filter(t => t.bounds.geoOverlaps(xmin, ymin, xmax, ymax))
    if (faces.isEmpty) throw new Exception("Triangle query did not match any octahedral faces!")

    // recurse, and combine results
    faces.tail.foldLeft(faces.head.getRangesRecursively(xmin, ymin, xmax, ymax, maxDepth))((acc, face) => {
      val faceRanges: Seq[IndexRange] = face.getRangesRecursively(xmin, ymin, xmax, ymax, maxDepth)
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
    val bounds: TriBounds = TriBounds(Extent(x0, x0 + 90.0), Extent(x0, x0 + 90.0), Extent(x0, x0 + 90.0), Extent(y0, y0 + 90.0))
    Triangle(idx, orientation, bounds, 1)
  }

  def geoToOctX(geoX: Double, y: Double): Double = {
    // the octahedral faces are all 90 degrees wide, so you just need to pick the
    // right 90-degree range that includes the given geoX
    val geoX0: Double = 90.0 * Math.floor((180.0 + Math.max(-180.0, Math.min(179.999999, geoX))) / 90.0) - 180.0

    // given the range minimum, transform the geo-X into an oct-X
    val pY = (90.0 - Math.abs(y)) / 90.0
    val xMid = geoX0 + 45.0
    val octX = xMid - pY * (xMid - geoX)
    octX
  }

  def octToGeoX(octX: Double, y: Double): Double = {
    // the octahedral faces are all 90 degrees wide, so you just need to pick the
    // right 90-degree range that includes the given geoX
    val geoX0: Double = 90.0 * Math.floor((180.0 + Math.max(-180.0, Math.min(179.999999, octX))) / 90.0) - 180.0

    // given the range minimum, transform the geo-X into an oct-X
    val pY = (90.0 - Math.abs(y)) / 90.0
    // mind the potential division-by-zero
    if (pY < 1e-6) return octX
    val xMid = geoX0 + 45.0
    val geoX = xMid - (xMid - octX) / pY
    geoX
  }

  // assume you've given a GEOGRAPHIC X-coordinate
  def getTriangle(geoX: Double, y: Double, maxDepth: Int): Triangle = {
    require(maxDepth > 0, s"maxDepth must be at least one; found $maxDepth")
    val octFace = getInitialTriangle(geoX, y)
    if (maxDepth > 1) {
      val result = octFace.getTriangle(geoToOctX(geoX, y), y, maxDepth - 1)
      result
    } else {
      val result = octFace
      result
    }
  }

  def indexBinaryString(index: Long, depth: Int): String =
    index.toBinaryString.reverse.padTo(3 * depth, "0").reverse.mkString("")

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
    val bounds: TriBounds = TriBounds(Extent(x0, x0 + 90.0), Extent(x0, x0 + 90.0), Extent(x0, x0 + 90.0), Extent(y0, y0 + 90.0))
    val octFace: Triangle = Triangle(octIdx, topOrientation(octXi.toInt, octYi.toInt), bounds, 1)

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
    val octX = geoToOctX(geoX, y)
    val x2 = octToGeoX(octX, y)
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

  // dump files out that we want to visualize
  // (these are not unit tests so much as manual validation tests)
  try {
    pw = new PrintWriter(new FileWriter("test-points.txt"))
    pw.println("name\tgeo_wkt\toct_wkt")
    LocationsByName.foreach {
      case (name, LatLon(y, geoX)) =>
        val t = getInitialTriangle(geoX, y)
        pw.println(s"$name\tPOINT($geoX $y)\tPOINT(${geoToOctX(geoX, y)} $y)")
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
      pw.println(s"${OrientationNames(t.orientation)}\t${t.index}\t${t.bitString}\t${t.octWkt}\t${t.geoWkt}\t${if (t.isAllApex) 1 else 0}")
    }
    pw.close()

    pw = new PrintWriter(new FileWriter("test-progression.txt"))
    pw.println("index_from\tindex_to\toct_wkt\tgeo_wkt")
    for (t <- iterator(4).sliding(2, 1)) t match {
      case Seq(t0, t1) =>
        pw.println(s"${t0.index}\t${t1.index}\t${geoLineWkt(t0, t1)}\t${octLineWkt(t0, t1)}")
    }
    pw.close()

  } finally {
    pw.close()
  }
}