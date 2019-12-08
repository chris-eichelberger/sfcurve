package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve.{GapMergedIndexRange, IndexRange, RangeComputeHints, SpaceFillingCurves}
import org.scalatest.{FunSpec, Matchers}

class ZCurve3DSpec extends FunSpec with Matchers {
  describe("Z3 as an SFC implementation") {
    it("Emits query ranges in a consistent order") {
      val sfc = new ZCurve3D(8)
      val hints = new RangeComputeHints
      hints.put(GapMergedIndexRange.HintsKeyMapGap, 1L.asInstanceOf[AnyRef])

      def validateQueryRanges(x0: Double, x1: Double, y0: Double, y1: Double, z0: Double, z1: Double): Boolean = {
        val ranges: Seq[IndexRange] = sfc.toRanges(x0, y0, z0, x1, y1, z1, Option(hints))
        println(f"Number of 3D Z-order ranges satisfying query X[$x0%1.4f, $x1%1.4f] Y[$y0%1.4f, $y1%1.4f] Z[$z0%1.4f, $z1%1.4f]:  ${ranges.size}%d")
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
      validateQueryRanges(-178.0, 179.0, -86.0, 87.0, 0.0, 1000.0) should be(true)

      // static query on a small area
      validateQueryRanges(-178.0, -178.0 + 1e-10, -86.0, 87.0, 0.0, 1000.0) should be(true)

      // generator queries to validate
      for (trial <- 1 to 100) {
        val x0: Double = -180.0 + 1e-10 + Math.random() * 350.0
        val x1: Double = x0 + 1.0 + (178.0 - x0) * Math.random()
        val y0: Double = -90.0 + 1e-10 + Math.random() * 170.0
        val y1: Double = y0 + 1.0 + (88.0 - y0) * Math.random()
        val z0: Double = sfc.zDimension.extent.min + 1e-10 + 0.9 * (sfc.zDimension.extent.max - sfc.zDimension.extent.min) * Math.random
        val z1: Double = z0 + 1.0 + 0.9 * (sfc.zDimension.extent.max - z0) * Math.random()
        (x0 < x1) should be(true)
        (y0 < y1) should be(true)
        (z0 < z1) should be(true)
        (x0 >= -180.0 && x0 < 180.0) should be(true)
        (x1 >= -180.0 && x1 < 180.0) should be(true)
        (y0 >= -90.0 && y0 < 90.0) should be(true)
        (y1 >= -90.0 && y1 < 90.0) should be(true)
        (z0 >= sfc.zDimension.extent.min && z0 < sfc.zDimension.extent.max) should be(true)
        (z1 >= sfc.zDimension.extent.min && z1 < sfc.zDimension.extent.max) should be(true)
        validateQueryRanges(x0, x1, y0, y1, z0, z1) should be(true)
      }
    }

  }

}
