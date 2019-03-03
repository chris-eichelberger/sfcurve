package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve.SpaceFillingCurves
import org.scalatest.{FunSpec, Matchers}

class ZCurve2DSpec extends FunSpec with Matchers {
  describe("SPI access") {
    it("provides an SFC") {
      val sfc = SpaceFillingCurves("zorder", Map(ZOrderSFCProvider.RESOLUTION_PARAM -> Int.box(100)))
      sfc.isInstanceOf[ZCurve2D] should equal(true)
    }

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
  }

}
