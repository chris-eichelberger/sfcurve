package org.locationtech.sfcurve

import java.text.SimpleDateFormat
import java.util.Date

import org.locationtech.sfcurve.Dimensions._
import org.locationtech.sfcurve.Utilities.CartesianProductIterable
import org.scalatest.{FunSpec, Matchers}

import scala.util.{Failure, Success, Try}

class DimensionsSpec extends FunSpec with Matchers {
  val dtf = new SimpleDateFormat("yyyy-MM-dd HH:mm")
  val dt0: Date = dtf.parse("1988-02-12 11:45")
  val dt1: Date = dtf.parse("1998-08-28 22:15")
  val dt2: Date = dtf.parse("2001-04-23 12:30")
  val dt3: Date = dtf.parse("2019-02-12 08:11")

  describe("date ordering") {
    it("should order dates properly") {
      dateOrdering.equiv(dt3, dt3) should equal(true)
      dateOrdering.equiv(dt2, dt3) should equal(false)
      dateOrdering.equiv(dt3, dt2) should equal(false)

      dateOrdering.lt(dt0, dt3) should equal(true)
      dateOrdering.lt(dt3, dt0) should equal(false)
      dateOrdering.lt(dt3, dt3) should equal(false)

      dateOrdering.gt(dt3, dt0) should equal(true)
      dateOrdering.gt(dt0, dt3) should equal(false)
      dateOrdering.gt(dt3, dt3) should equal(false)

      dateOrdering.gteq(dt3, dt0) should equal(true)
      dateOrdering.gteq(dt3, dt3) should equal(true)
      dateOrdering.gteq(dt0, dt3) should equal(false)

      dateOrdering.lteq(dt0, dt3) should equal(true)
      dateOrdering.lteq(dt3, dt3) should equal(true)
      dateOrdering.lteq(dt3, dt0) should equal(false)
    }
  }

  describe("Extent") {
    it("should instantiate different types") {
      Extent[Long](0L, 63L) should not be null
      Extent[Int](0, 63) should not be null
      Extent[Double](0.1, 63.7) should not be null
      Extent[Date](dt0, dt3) should not be null
    }

    it("should bound elements properly") {
      val outerDate = Extent(dt0, dt2, incMin = true, incMax = false)
      outerDate.contains(dt0) should equal(true)
      outerDate.contains(dt1) should equal(true)
      outerDate.contains(dt2) should equal(false)
      outerDate.contains(dt3) should equal(false)

      val outerLong = Extent(0L, 100L)
      outerLong.contains(-1L) should equal(false)
      (0 until 100).forall(i => outerLong.contains(i.toLong)) should equal(true)
      outerLong.contains(100L) should equal(false)
    }
  }

  describe("standard geographic dimensions") {
    val InvalidBin = -999

    it("should work as expected at various precisions") {
      Seq(2, 3, 11, 1024).foreach(cardinality => {
        println(s"Testing geographic dimensions with $cardinality bins...")
        val dimLongitude = Longitude(cardinality)
        val dimLatitude= Latitude(cardinality)

        def getBin(dim: Dimension[Double], value: Double): Long = {
          val bin = Try { dim.index(Seq(value)) } match {
            case Success(b) => b
            case Failure(_) => InvalidBin
          }
          println(s"  ${dim.extent} given $value -> bin $bin over 0..${cardinality - 1}")
          bin
        }

        dimLongitude.cardinality should equal(cardinality)
        getBin(dimLongitude, -540.0) should equal(0)
        getBin(dimLongitude, -180.0) should equal(0)
        getBin(dimLongitude, 0.0) should equal(cardinality  >> 1)
        getBin(dimLongitude, 180.0 - 1e-10) should equal(cardinality - 1)
        getBin(dimLongitude, 180.0) should equal(0)  // Longitude DOES wrap at 180!
        getBin(dimLongitude, 360.0) should equal(cardinality >> 1)

        dimLatitude.cardinality should equal(cardinality)
        getBin(dimLatitude, -270.0) should equal(0)
        getBin(dimLatitude, -90.0) should equal(0)
        getBin(dimLatitude, 0.0) should equal(cardinality >> 1)
        getBin(dimLatitude, 90.0 - 1e-10) should equal(cardinality - 1)
        getBin(dimLatitude, 90.0) should equal(cardinality - 1)  // Latitude does NOT wrap at 90!
        getBin(dimLatitude, 180.0) should equal(cardinality >> 1)
      })
    }
  }

  describe("two-dimensional cell") {
    it("should compute proper containment") {
      val lonExtent: Extent[_] = Longitude(10).extent
      val latExtent: Extent[_] = Latitude(10).extent
      val cell2d = Cell(Vector(lonExtent, latExtent))

      cell2d.contains(Vector(-180.0, -90.0)) should equal(true)
      cell2d.contains(Vector(0.0, 0.0)) should equal(true)
      cell2d.contains(Vector(180.0, 0.0)) should equal(false)
      cell2d.contains(Vector(0.0, 90.0)) should equal(true)
      cell2d.contains(Vector(180.0, 90.0)) should equal(false)
    }

    it("should add properly") {
      val exA = Extent[Long](10L, 19L)
      val exB = Extent[Double](-13.3, -9.7)

      val cellA = Cell(Vector(exA))
      cellA.contains(Vector(9L)) should equal(false)
      cellA.contains(Vector(10L)) should equal(true)
      cellA.contains(Vector(19L)) should equal(false)
      cellA.contains(Vector(20L)) should equal(false)

      val cellB = Cell(Vector(exB))
      cellB.contains(Vector(-13.4)) should equal(false)
      cellB.contains(Vector(-13.3)) should equal(true)
      cellB.contains(Vector(-9.8)) should equal(true)
      cellB.contains(Vector(-9.7)) should equal(false)
      cellB.contains(Vector(1.0)) should equal(false)

      val cellC = cellA + cellB

      // observe that the `Vector[Any]` annotation is required for this to pass;
      // without it, plain `Vector`, Scala converts all elements to type Double
      cellC.contains(Vector[Any](9L, -13.4)) should equal(false)
      cellC.contains(Vector[Any](10L, -13.4)) should equal(false)
      cellC.contains(Vector[Any](9L, -13.3)) should equal(false)
      cellC.contains(Vector[Any](10L, -13.3)) should equal(true)
    }
  }

  describe("generic space-filling curve") {
    // dummy class for testing
    def R(subs: Vector[Discretizor], maxGap: Long = 0): RowMajorSFC = RowMajorSFC(subs, maxGap)

    it("should work for simple 2D row-major form on a WGS84 plane") {
      val r = R(Vector(Longitude(8), Latitude(4)))

      def test(x: Double, y: Double, expectedIndex: Long, xi: Long, yi: Long): Unit = {
        val actualIndex: Long = r.index(Seq(x, y))
        println(s"R(8,4) at POINT($x $y))...")
        println(s"  fold:  ($xi, $yi) -> $actualIndex (expected $expectedIndex)")
        actualIndex should equal(expectedIndex)
        val cell = r.inverseIndex(expectedIndex)
        val coord = r.unfold(expectedIndex)
        println(s"  unfold:  $expectedIndex -> (${coord(0)}, ${coord(1)})")
        println(s"  unfold:  $expectedIndex -> $cell")
        coord(0) should equal(xi)
        coord(1) should equal(yi)
        cell.contains(Vector(x, y)) should equal(true)
      }

      // test fold and unfold
      test(-180.0, -90.0, 0L, 0, 0)
      test(179.9, -90.0, 28, 7, 0)
      test(0.0, 0.0, 18, 4, 2)
      test(-180.0, 89.9, 3, 0, 3)
      test(179.9, 89.9, 31L, 7, 3)

      // test query ranges...
      val hints = new RangeComputeHints()
      hints.put("MaximumAllowableGapBetweenQueryRanges", 0L.asInstanceOf[AnyRef])
      val extX = Extent[Double](0.0, 180.0)
      val extY = Extent[Double](0.0, 90.0)

      val h2x: Seq[IndexRange] = r.queryRanges(Seq(Option(extX), None), Option(hints)).toList
      println(s"R(8,4) query-ranges($extX, None):  ${h2x.map(ir => if (ir.lower == ir.upper) ir.lower.toString else s"${ir.lower}-${ir.upper}").mkString("[", ", ", "]")}")
      h2x.size should equal(1)
      h2x.head.lower should equal(16)
      h2x.head.upper should equal(31)

      val h2xh2y: Seq[IndexRange] = r.queryRanges(Seq(Option(extX), Option(extY)), Option(hints)).toList
      println(s"R(8,4) query-ranges($extX, $extY):  ${h2xh2y.map(ir => if (ir.lower == ir.upper) ir.lower.toString else s"${ir.lower}-${ir.upper}").mkString("[", ", ", "]")}")
      h2xh2y.size should equal(4)
      h2xh2y.head.lower should equal(18)
      h2xh2y.head.upper should equal(19)

      hints.put("MaximumAllowableGapBetweenQueryRanges", 5L.asInstanceOf[AnyRef])
      val h2xh2ym: Seq[IndexRange] = r.queryRanges(Seq(Option(extX), Option(extY)), Option(hints)).toList
      println(s"R(8,4) query-ranges($extX, $extY):  ${h2xh2ym.map(ir => if (ir.lower == ir.upper) ir.lower.toString else s"${ir.lower}-${ir.upper}").mkString("[", ", ", "]")}")
      h2xh2ym.size should equal(1)
      h2xh2ym.head.lower should equal(18)
      h2xh2ym.head.upper should equal(31)
    }
  }
}
