package org.locationtech.sfcurve.zorder

import java.io.{FileWriter, PrintWriter}

import org.locationtech.sfcurve.Dimensions.{Cell, Extent, Latitude, Longitude}
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

case class Triangle(index: Long, orientation: Int, X: Extent[Double], Y: Extent[Double], depth: Int) {
  val x0: Double = X.min
  val x1: Double = X.max
  val xMid: Double = 0.5 * (x0 + x1)
  val x14: Double = 0.75 * x0 + 0.25 * x1
  val x34: Double = 0.25 * x0 + 0.75 * x1
  val (y0: Double, y1: Double, mLeftInv: Double, mRightInv: Double) = if (isApexUp(orientation)) {
    (Y.min, Y.max, MNegInv, MPosInv)
  } else {
    (Y.max, Y.min, MPosInv, MNegInv)
  }
  val yMid: Double = 0.5 * (y0 + y1)

  def bitString: String = indexBinaryString(index, depth)

  def stackIndex(nextIndex: Long): Long = (index << 3) | (nextIndex & 7)

  def nextOrientation(transition: Int): Int = OrientationTransitions((orientation, transition))

  def child(transition: Int): Triangle = transition match {
    case TransCenter => Triangle(stackIndex(TransCenter), nextOrientation(TransCenter), Extent(x14, x34), makeExtent(y0, yMid), depth + 1)
    case TransApex =>   Triangle(stackIndex(TransApex), nextOrientation(TransApex), Extent(x14, x34), makeExtent(yMid, y1), depth + 1)
    case TransLL =>     Triangle(stackIndex(TransLL), nextOrientation(TransLL), Extent(x0, xMid), makeExtent(y0, yMid), depth + 1)
    case TransLR =>     Triangle(stackIndex(TransLR), nextOrientation(TransLR), Extent(xMid, x1), makeExtent(y0, yMid), depth + 1)
    case _ =>           throw new Exception(s"Invalid transition $transition")
  }

  def getTriangle(x: Double, y: Double, maxDepth: Int): Triangle = {
    require(X.contains(x), s"X is out of bounds:  $x notIn $X")
    require(Y.contains(y), s"Y is out of bounds:  $y notIn $Y")

    if (maxDepth < 1) return this

    // apex
    if (Math.abs(y - y1) < Math.abs(y - y0)) {
      require(x >= x14 && x <= x34, s"Cannot recurse into apex, X $x is out of range ($x14, $x34)")
      return child(TransApex).getTriangle(x, y, maxDepth - 1)
    }

    // you know you're recursing into the base

    // check left and right
    if (x < xMid) {
      // check left
      if (x < x14) {
        return child(TransLL).getTriangle(x, y, maxDepth - 1)
      }
      val xProbe: Double = mLeftInv * (y - y0) + xMid
      if (xProbe >= x0 && xProbe < xMid && x < xProbe) {
        return child(TransLL).getTriangle(x, y, maxDepth - 1)
      }
    } else {
      // check right
      if (x > x34) {
        return child(TransLR).getTriangle(x, y, maxDepth - 1)
      }
      val xProbe: Double = mRightInv * (y - y0) + xMid
      if (xProbe >= xMid && xProbe <= x1 && x >= xProbe) {
        return child(TransLR).getTriangle(x, y, maxDepth - 1)
      }
    }

    // if you get this far, recurse center
    child(TransCenter).getTriangle(x, y, maxDepth - 1)
  }

  def polygonPoints: String = orientation match {
    case A_BC | A_CB | B_AC | B_CA | C_AB | C_BA =>
      // apex on top
      s"${xMid} ${Y.max}, ${X.max} ${Y.min}, ${X.min} ${Y.min}, ${xMid} ${Y.max}"
    case AB_C | AC_B | BA_C | BC_A | CA_B | CB_A =>
      // apex on bottom
      s"${X.min} ${Y.max}, ${X.max} ${Y.max}, ${xMid} ${Y.min}, ${X.min} ${Y.max}"
    case _ =>
      throw new Exception(s"Invalid orientation $orientation")
  }

  def wkt: String = {
    "POLYGON((" + polygonPoints + "))"
  }
}

object TriN {
  val TriTopLongitude = Longitude(4)
  val TriTopLatitude = Latitude(2)

  // possible triangle orientations...
  // apex up
  val A_BC = 0
  val B_CA = 1
  val C_AB = 2
  val A_CB = 3
  val B_AC = 4
  val C_BA = 5
  // apex down
  val AB_C = 6
  val BC_A = 7
  val CA_B = 8
  val AC_B = 9
  val BA_C = 10
  val CB_A = 11

  val Orientations: Set[Int] = Set(
    A_BC, B_CA, C_AB, A_CB, B_AC, C_BA,
    AB_C, BC_A, CA_B, AC_B, BA_C, CB_A
  )

  val TransCenter = 0
  val TransApex = 1
  val TransLR = 2
  val TransLL = 4

  val Transitions: Set[Int] = Set(TransCenter, TransApex, TransLR, TransLL)

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

  def topOrientation(ix: Int, iy: Int): Int = {
    require(ix >= 0 && ix < 4)
    require(iy >= 0 && iy < 2)
    (iy & 1, ix & 1) match {
      case (0, 0) => 11
      case (0, 1) => 7
      case (1, 0) => 0
      case (1, 1) => 3
    }
  }

  // map the given (x, y) coordinate to one face of the octahedron
  // that wraps the entire (spherical) Earth
  def getInitialTriangle(x: Double, y: Double): Triangle = {
    val ix = TriTopLongitude.index(Seq(x))
    val iy = TriTopLatitude.index(Seq(y))
    val orientation: Int = topOrientation(ix.toInt, iy.toInt)
    val idx: Long = (iy << 2) | (if (ix < 2) ix else 5 - ix)
    val x0: Double = -180.0 + 90.0 * ix
    val y0: Double = -90.0 + 90.0 * iy
    Triangle(idx, orientation, Extent(x0, x0 + 90.0), Extent(y0, y0 + 90.0), 1)
  }

  def getXT(x: Double, y: Double, X: Extent[Double]): Double = {
    val dx = X.max - X.min
    val pY = (90.0 - Math.abs(y)) / 90.0
    val xMid = X.min + 0.5 * dx
    val xT = xMid - pY * (xMid - x)
    xT
  }

  def getXTInv(xT: Double, y: Double, X: Extent[Double]): Double = {
    val dx = X.max - X.min
    val pY = (90.0 - Math.abs(y)) / 90.0
    val xMid = X.min + 0.5 * dx
    val x = if (pY > 1e-6) xMid - (xMid - xT) / pY else xT
    x
  }

  def getTriangle(x: Double, y: Double, maxDepth: Int): Triangle = {
    require(maxDepth > 0, s"maxDepth must be at least one; found $maxDepth")
    val t0 = getInitialTriangle(x, y)
    if (maxDepth > 1) {
      val result = t0.getTriangle(getXT(x, y, t0.X), y, maxDepth - 1)
      result
    } else {
      val result = t0
      result
    }
  }

  def indexBinaryString(index: Long, depth: Int): String =
    index.toBinaryString.reverse.padTo(3 * depth, "0").reverse.mkString("")

  def index(x: Double, y: Double, maxDepth: Int): Long = {
    getTriangle(x, y, maxDepth).index
  }

  def invIndex(idx: Long, depth: Int): Triangle = {
    require(depth > 0, s"Depth ($depth) must be at least 1")

    val octIdx = (idx >> (3 * (depth - 1))) & 7
    val octIdxX = octIdx & 3
    val octXi = if (octIdxX < 2) octIdxX else 5 - octIdxX
    val octYi = (octIdx >> 2) & 1
    val x0 = octXi * 90.0 - 180.0
    val y0 = octYi * 90.0 - 90.0
    val t0: Triangle = Triangle(octIdx, topOrientation(octXi.toInt, octYi.toInt), Extent(x0, x0 + 90.0), Extent(y0, y0 + 90.0), 1)

    println(s"invIndex(${indexBinaryString(idx, depth)}, $depth)...")

    (depth - 2  to 0 by -1).foldLeft(t0)((acc, i) => {
      val triad = (idx >> (3 * i)) & 7
      println(s"  triad $i:  ${indexBinaryString(triad, 1)}")
      acc.child(triad.toInt)
    })
  }

  def iterator(depth: Int): Iterator[Triangle] = {
    require(depth > 0, s"Depth ($depth) must be a positive integer")

    val faces: Seq[Int] = (0 to 7).toSeq
    val seqs: Seq[Seq[Int]] = (2 to depth).foldLeft(Seq(faces))((acc, _) => acc ++ Seq(Transitions.toSeq))

    CartesianProductIterable(seqs).iterator.map( seq => {
      val index = seq.foldLeft(0L)((acc, i) => i match {
        case n: Int =>
          (acc << 3) | n
      })
      val t = invIndex(index, depth)
      println(s"  iteration:  ${seq.mkString("(", ", ", ")")} -> ${indexBinaryString(index, depth)} == ${t.bitString}")
      require(t.index == index)
      t
    })
  }
}

object NamedLocations {
  case class LatLon(latitude: Double, longitude: Double)
  val LocationsByName: Map[String, LatLon] = Map(
    "Charlottesville" -> LatLon(38.0293, -78.4767),
    "Eichelberger" -> LatLon(38.054444, -78.688256),
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
    "Gisborne" -> LatLon(-38.670242, 178.026712),
    "South Pole" -> LatLon(-90.0, -179.0),
    "North Pole" -> LatLon(90.0 - 1e-10, 179.0),
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

  def testInitialTri(n: Int, x: Double, y: Double): Unit = {
    val t = getTriangle(x, y, 1)
    val bIndex = t.index.toBinaryString.reverse.padTo(3, "0").reverse.mkString("")
    val xT = getXT(x, y, t.X)
    val x2 = getXTInv(xT, y, t.X)
    pw.println(f"$n%d\tPOINT($x%1.4f $y%1.4f)\tPOINT($xT%1.4f $y%1.4f)\t$bIndex%s")
    if (Math.abs(y) < 89.5) {
      require(Math.abs(x2 - x) <= 1e-6, f"Failed to satisfy XT inverse:  $x%1.6f <> $xT%1.6f")
    }
  }

  try {
    pw = new PrintWriter(new FileWriter("test-points.txt"))
    pw.println("name\torig_wkt\toct_wkt")
    LocationsByName.foreach {
      case (name, LatLon(y, x)) =>
        val t = getInitialTriangle(x, y)
        pw.println(s"$name\tPOINT($x $y)\tPOINT(${getXT(x, y, t.X)} $y)")
    }
    pw.close()

    pw = new PrintWriter(new FileWriter("test-triangles.txt"))
    pw.println("index\twkt")
    for (row <- 0 to 1; col <- 0 to 3) {
      val x0 = -180.0 + 90.0 * col
      val y0 = -90.0 + 90.0 * row
      val x1 = x0 + 90.0 - 1e-6
      val y1 = y0 + 90.0 - 1e-6
      val xMid = 0.5 * (x0 + x1)
      val yMid = 0.5 * (x0 + x1)
      val index = getInitialTriangle(xMid, yMid).index.toBinaryString.reverse.padTo(3, "0").reverse.mkString("")
      if (row == 0) {
        pw.println(s"$index\tPOLYGON(($x0 $y1, $xMid $y0, $x1 $y1, $x0 $y1))")
      } else {
        pw.println(s"$index\tPOLYGON(($x0 $y0, $xMid $y1, $x1 $y0, $x0 $y0))")
      }
    }
    pw.close()

    pw = new PrintWriter(new FileWriter("test.txt"))
    pw.println("label\torig_wkt\ttri_wkt\tindex")
    val xs = (-180.0 to 179.999 by 22.5).zipWithIndex
    val ys = (-90.0 to 89.999 by 11.25*0.5).zipWithIndex
    val ny = ys.length
    for (xx <- xs; yy <- ys) {
      testInitialTri(yy._2 * ny + xx._2, xx._1, yy._1)
    }
    pw.close()

    pw = new PrintWriter(new FileWriter("test-index.txt"))
    pw.println("depth\tindex_dec\tindex_bits\twkt")
    val target: LatLon = LocationsByName("Eichelberger")
    (1 to 21).foldLeft((target.longitude, target.latitude))((acc, depth) => acc match {
      case (x, y) =>
        val t = getTriangle(target.longitude, target.latitude, depth)
        pw.println(depth + "\t" + t.index + "\t" + t.bitString + "\t" + t.wkt)
        (x, y)
    })
    pw.close()

    pw = new PrintWriter(new FileWriter("test-tiles.txt"))
    pw.println("index_dec\tindex_bits\twkt")
    for (t <- iterator(2)) {
      pw.println(s"${t.index}\t${t.bitString}\t${t.wkt}")
    }
  } finally {
    pw.close()
  }
}