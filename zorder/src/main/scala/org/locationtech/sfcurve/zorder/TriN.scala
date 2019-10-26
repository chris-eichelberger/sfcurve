package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve.Dimensions.{Cell, Extent, Latitude, Longitude}

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

case class Triangle(index: Long, orientation: Int, X: Extent[Double], Y: Extent[Double]) {
  def getTriangle(x: Double, y: Double, maxDepth: Int): Triangle = {
    if (maxDepth < 1) return this
    
    val x0 = X.min
    val x1 = X.max
    val dx = x1 - x0
    val x14 = x0 + 0.25 * dx
    val xmid = x0 + 0.50 * dx
    val x34 = x0 + 0.75 * dx

    val y0 = Y.min
    val y1 = Y.max
    val dy = y1 - y0
    val ymid = y0 + 0.50 * dy

    val iy: Int = if (y >= 0.5 * (Y.min + Y.max)) 1 else 0
      
    // two cases:  apex up; or apex down

    // TODO!
    null
  }
}

object TriN {
  val Sqrt3: Double = Math.sqrt(3.0)
  val TriTopLongitude = Longitude(4)
  val TriTopLatitude = Latitude(2)

  // map the given (x, y) coordinate to one face of the octahedron
  // that wraps the entire (spherical) Earth
  def getInitialTriangle(x: Double, y: Double): Triangle = {
    val ix = TriTopLongitude.index(Seq(x))
    val iy = TriTopLatitude.index(Seq(y))
    val orientation: Int = (iy & 1, ix & 1) match {
      case (0, 0) => 11
      case (0, 1) => 7
      case (1, 0) => 0
      case (1, 1) => 3
    }
    val idx: Long = (iy << 2) | (if (ix < 2) ix else 5 - ix)
    val x0: Double = -180.0 + 90.0 * ix
    val y0: Double = -90.0 + 90.0 * iy
    Triangle(idx, orientation, Extent(x0, x0 + 90.0), Extent(y0, y0 + 90.0))
  }

  def getXT(x: Double, y: Double, X: Extent[Double]): Double = {
    val dx = X.max - X.min
    val dxT = Math.cos(Math.toRadians(y))
    val pX0 = x / dx
    (dx - 0.5 * dxT) + dxT * pX0
  }

  def getTriangle(x: Double, y: Double, maxDepth: Int): Triangle = {
    require(maxDepth > 0, s"maxDepth must be at least one; found $maxDepth")
    val t0 = getInitialTriangle(x, y)
    if (maxDepth > 1) {
      t0.getTriangle(getXT(x, y, t0.X), y, maxDepth - 1)
    } else {
      t0
    }
  }

  def index(x: Double, y: Double, maxDepth: Int): Long = {
    getTriangle(x, y, maxDepth).index
  }

  def invIndex(idx: Long): Triangle = {
    // TODO!
    null
  }
}

class TriN {
}

object TriTest extends App {
  def testInitialTri(x: Double, y: Double): Unit = {
    val t = getTriangle(x, y, 1)
    val bIndex = t.index.toBinaryString.reverse.padTo(3, "0").reverse.mkString("")
    println(f"POINT($x%1.4f, $y%1.4f) -> T $bIndex%s")
  }

  testInitialTri(45.0, 45.0)
  testInitialTri(0.0, 0.0)
  testInitialTri(45.0, -45.0)
}