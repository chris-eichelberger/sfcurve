/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve.hilbert

import org.locationtech.sfcurve.Dimensions.{Latitude, Longitude}
import org.locationtech.sfcurve.{GapMergedIndexRange, IndexRange, RangeComputeHints, SpaceFillingCurves}
import org.scalatest.MustMatchers.convertToAnyMustWrapper
import org.scalatest._

class HilbertCurveSpec extends FunSpec with Matchers {

  val EPSILON: Double = 1E-3

  describe("SPI access") {
    it("resolves a HilbertCurve2D") {
      val sfc = SpaceFillingCurves("hilbert", Map(HilbertCurve2DProvider.RESOLUTION_PARAM -> Int.box(10)))
      sfc.isInstanceOf[HilbertCurve2D] should equal(true)
    }
  }
  describe("A HilbertCurve implementation using UG lib") {

    it("translates (Double,Double) to Long and Long to (Double, Double)"){
      val bitsPerDim: Int = 16  // 32 bits total in the curve
      val gridCellsPerDim: Long = 1L << bitsPerDim

      for (atX <- Seq(-180.0, -179.9999, -10, 0, 10, 179.9999, 180); atY <- Seq(-90, -89.9999, -10, 0, 10, 89.9999, 90)) {
        val sfc = new HilbertCurve2D(bitsPerDim)  // 32 bits (4 giga cells) in the curve
        val index: Long = sfc.toIndex(atX, atY)

        val point = sfc.toPoint(index)

        val dimLongitude = Longitude(gridCellsPerDim.toLong)
        val dimLatitude = Latitude(gridCellsPerDim.toLong)

        dimLongitude.toExtent(dimLongitude.toBin(atX)).contains(point._1) should be(true)
        dimLatitude.toExtent(dimLatitude.toBin(atY)).contains(point._2) should be(true)
      }
    }

    it("implements a range query"){

      val sfc = new HilbertCurve2D(3)  // 6 bits (64 cells) in the curve

      // any gaps of only one index cell will be merged automatically
      val hints = new RangeComputeHints
      hints.put(GapMergedIndexRange.HintsKeyMapGap, 1L.asInstanceOf[AnyRef])

      val range = sfc.toRanges(-178.123456, -86.398493, 179.3211113, 87.393483, Option(hints))

      range should have length 6

      // the last range is not wholly contained within the query region
      val (_, _, lastcontains) = range(2).tuple
      lastcontains should be(false)
    }

    it("Takes a Long value to a Point (Double, Double)"){

      val value = 0L
      val sfc = new HilbertCurve2D(8)  // 16 bits (65536 cells) in the curve
      val point: (Double, Double) = sfc.toPoint(value)
      print(point)     

    }

    it("Takes a Point (Double, Double) to a Long value"){

      val sfc = new HilbertCurve2D(8)  // 16 bits (65536 cells) in the curve
      val value: Long = sfc.toIndex(0.0,0.0)
      print(value)     

    }

    it("Emits query ranges in a consistent order") {
      val sfc = new HilbertCurve2D(8)  // 16 bits (65536 cells) in the curve
      val hints = new RangeComputeHints
      hints.put(GapMergedIndexRange.HintsKeyMapGap, 1L.asInstanceOf[AnyRef])

      def validateQueryRanges(x0: Double, x1: Double, y0: Double, y1: Double): Boolean = {
        val ranges: Seq[IndexRange] = sfc.toRanges(x0, y0, x1, y1, Option(hints))
        println(f"Number of Hilbert ranges satisfying query X[$x0%1.4f, $x1%1.4f] Y[$y0%1.4f, $y1%1.4f]:  ${ranges.size}%d")
        ranges.sliding(2, 1).foreach {
          case Seq(rangePrevious: IndexRange, rangeNext: IndexRange) =>
            if (rangePrevious.upper <= rangeNext.lower) {
              println(s"Sliding Hilbert range pair FAILED:  [${rangePrevious.lower}, ${rangePrevious.upper}], [${rangeNext.lower}, ${rangeNext.upper}]")
              return false
            }
        }
        true
      }

      // static query to validate (355 ranges)
      validateQueryRanges(-178.0, 179.0, -86.0, 87.0) must be(true)

      // generator queries to validate
      for (trial <- 1 to 100) {
        // NB:  Uzaygezen does NOT like it when you query on too small a dimension difference
        // (it gets unhappy when the X or Y indexes are not different, i.e., [17, 17] is bad)
        val x0: Double = -180.0 + 1e-10 + Math.random() * 350.0
        val x1: Double = x0 + 1.0 + (178.0 - x0) * Math.random()
        val y0: Double = -90.0 + 1e-10 + Math.random() * 170.0
        val y1: Double = y0 + 1.0 + (88.0 - y0) * Math.random()
        (x0 < x1) should be(true)
        (y0 < y1) should be(true)
        (x0 >= -180.0 && x0 < 180.0) should be(true)
        (x1 >= -180.0 && x1 < 180.0) should be(true)
        (y0 >= -90.0 && y0 < 90.0) should be(true)
        (y1 >= -90.0 && y1 < 90.0) should be(true)
        validateQueryRanges(x0, x1, y0, y1)
      }
    }

  }
}
