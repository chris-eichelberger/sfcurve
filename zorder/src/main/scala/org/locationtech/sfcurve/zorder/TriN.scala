package org.locationtech.sfcurve.zorder

import java.io.{FileWriter, PrintWriter}

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
    val pY = (90.0 - Math.abs(y)) / 90.0
    val xMid = X.min + 0.5 * dx
    val xT = xMid - pY * (xMid - x)
//    println(f"getXT($x%1.4f, $y%1.4f, [${X.min}%1.4f, ${X.max}%1.4f])...")
//    println(f"  dx $dx%1.4f")
//    println(f"  xMid $xMid%1.4f")
//    println(f"  pY $pY%1.4f")
//    println(f"  xMid - xT")
//    println(f"  xT $xT%1.4f")
    xT
  }

  def getXTInv(xT: Double, y: Double, X: Extent[Double]): Double = {
    val dx = X.max - X.min
    val pY = (90.0 - Math.abs(y)) / 90.0
    val xMid = X.min + 0.5 * dx
    val x = if (pY > 1e-6) xMid - (xMid - xT) / pY else xT
//    println(f"getXTInv($xT%1.4f, $y%1.4f, [${X.min}%1.4f, ${X.max}%1.4f])...")
//    println(f"  dx $dx%1.4f")
//    println(f"  xMid $xMid%1.4f")
//    println(f"  pY $pY%1.4f")
//    println(f"  xMid - xT")
//    println(f"  x $x%1.4f")
    x
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
  var pw: PrintWriter = null

  def testInitialTri(n: Int, x: Double, y: Double): Unit = {
    val t = getTriangle(x, y, 1)
    val bIndex = t.index.toBinaryString.reverse.padTo(3, "0").reverse.mkString("")
    val xT = getXT(x, y, t.X)
//    println(f"  x' $xT%1.4f")
    val x2 = getXTInv(xT, y, t.X)
//    println(f"  x2 $x2%1.4f")
    pw.println(f"$n%d\tPOINT($x%1.4f $y%1.4f)\tPOINT($xT%1.4f $y%1.4f)\t$bIndex%s")
    require(Math.abs(x2 - x) <= 1e-6, f"Failed to satisfy XT inverse:  $x%1.6f <> $xT%1.6f")
  }

  try {
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
    val xs = (-179.999 to 179.999 by 22.5).zipWithIndex
    val ys = (-89.999 to 89.999 by 11.25*0.5).zipWithIndex
    val ny = ys.length
    for (xx <- xs; yy <- ys) {
      testInitialTri(yy._2 * ny + xx._2, xx._1, yy._1)
    }
//    testInitialTri(1, 45.0, 45.0)
//    testInitialTri(2, 0.0, 0.0)
//    testInitialTri(3, 45.0, -45.0)
//    testInitialTri(4, -180.0, -45.0)
//    testInitialTri(5, 179.9999, 45.0)
//    testInitialTri(6, 0.0, 89.9999)
//    testInitialTri(7, -78.4767, 38.0293)
  } finally {
    pw.close()
  }
}