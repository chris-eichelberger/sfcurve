package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve.{GapMergedIndexRange, IndexRange, RangeComputeHints, SpaceFillingCurves}
import org.scalatest.{FunSpec, Matchers}

class ZCurve2DSpec extends FunSpec with Matchers {
  describe("SPI access") {
    it("provides an SFC") {
      val sfc = SpaceFillingCurves("zorder", Map(ZOrderSFCProvider.RESOLUTION_PARAM -> Int.box(100)))
      sfc.isInstanceOf[ZCurve2D] should equal(true)
    }
  }

  describe("Z2 as an SFC implementation") {
    it("computes a covering set of ranges") {
      val sfc = SpaceFillingCurves("zorder", Map(ZOrderSFCProvider.RESOLUTION_PARAM -> Int.box(1024)))
      val ranges = sfc.toRanges(-80.0, 35.0, -75.0, 40.0, ZCurve2D.hints(maxRecurse = 32))

      ranges.length shouldBe 44
      println(s"Number of ZCurve2DSpec ranges:  ${ranges.length}")
      println(s"ZCurve2DSpec ranges:  ${ranges.mkString("[\n  ", ",\n  ", "\n]")}")
      val (l, r, contains) = ranges.head.tuple

      // This is the (top-quadrant) sketch of this unit test,
      // assuming a standard XYXY...XYXY bit interleaving:
      //
      //  180 +-----------+-----------+
      //      | 524288 to | 786432 to |
      //      | 786431    | 1048575   |
      //    0 +-----------+-----------+
      //      | 0 to      | 262144 to |
      //      | 262143    | 524287    |
      // -180 +-----------+-----------+
      //      -90         0           90

      // The test range (X[-80, -75], Y[35, 40]) should fall in the third quadrant
      // visited, falling in the index-range of 524288-786431.
      l shouldBe 631162
      r shouldBe 631163

      // old results:
      //l shouldBe 197616
      //r shouldBe 197631

      contains shouldBe true
    }

    it("Emits query ranges in a consistent order") {
      val sfc = new ZCurve2D(8)
      val hints = new RangeComputeHints
      hints.put(GapMergedIndexRange.HintsKeyMapGap, 1L.asInstanceOf[AnyRef])

      def validateQueryRanges(x0: Double, x1: Double, y0: Double, y1: Double): Boolean = {
        val ranges: Seq[IndexRange] = sfc.toRanges(x0, y0, x1, y1, Option(hints))
        println(f"Number of 2D Z-order ranges satisfying query X[$x0%1.4f, $x1%1.4f] Y[$y0%1.4f, $y1%1.4f]:  ${ranges.size}%d")
        val windows: List[Seq[IndexRange]] = ranges.sliding(2, 1).toList
        if (windows.size == 1) return true
        windows.foreach {
          case Seq(rangePrevious: IndexRange, rangeNext: IndexRange) =>
            if (rangePrevious.upper >= rangeNext.lower) {
              println(s"Sliding Z-order range pair FAILED:  [${rangePrevious.lower}, ${rangePrevious.upper}], [${rangeNext.lower}, ${rangeNext.upper}]")
              return false
            }
        }
        true
      }

      // static query to validate a large area
      validateQueryRanges(-178.0, 179.0, -86.0, 87.0) should be(true)

      // static query on a small area
      validateQueryRanges(-178.0, -178.0 + 1e-10, -86.0, 87) should be(true)

      // generator queries to validate
      for (trial <- 1 to 100) {
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
        validateQueryRanges(x0, x1, y0, y1) should be(true)
      }
    }

  }

}
