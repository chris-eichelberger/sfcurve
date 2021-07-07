package org.locationtech.sfcurve

import org.locationtech.sfcurve.Dimensions._
import org.locationtech.sfcurve.ImplicitCasts._
import org.scalatest.MustMatchers.convertToAnyMustWrapper
import org.scalatest.{FunSpec, Matchers}

import java.text.SimpleDateFormat
import java.util.Date

class ComposedCurveSpec extends FunSpec with Matchers {
  val dtf = new SimpleDateFormat("yyyy-MM-dd HH:mm")
  val dt0: Date = dtf.parse("1988-02-12 11:45")
  val dt1: Date = dtf.parse("1998-08-28 22:15")
  val dt2: Date = dtf.parse("2001-04-23 12:30")
  val dt3: Date = dtf.parse("2019-02-12 08:11")

  val BitsPerDimension: Int = 5
  val CardinalityPerDimension: Long = 1L << BitsPerDimension.toLong
  val X: Longitude = Longitude(CardinalityPerDimension)
  val Y: Latitude = Latitude(CardinalityPerDimension)
  val T: DateDimension = ExampleEra(CardinalityPerDimension)

  val qTOpt: Option[Extent[_]] = Option(Extent(dt0, dt3))
  val qXOpt: Option[Extent[_]] = Option(Extent(-79.0, -78.0))
  val qYOpt: Option[Extent[_]] = Option(Extent(38.0, 39.0))

  def R(children: Discretizor*): RowMajorSFC = RowMajorSFC(Vector(children:_*), 0)

  describe("composed-curve variants using only R, x, y, and t") {
    val Rtxy: SpaceFillingCurve = R(T, X, Y)
    val Rxy_t: SpaceFillingCurve = R(R(X, Y), T)
    val Rt_xy: SpaceFillingCurve = R(T, R(X, Y))
    val Rxt_y: SpaceFillingCurve = R(R(X, T), Y)

    it("must have been able to instantiate the curves") {
      Rt_xy must not be(null)
      Rxy_t must not be(null)
      Rtxy must not be(null)
      Rxt_y must not be(null)
    }

    it("must produce valid query ranges for different curves") {
      val queryTXY: Seq[Option[Extent[_]]] = Seq(qTOpt, qXOpt, qYOpt)
      val queryXYT: Seq[Option[Extent[_]]] = Seq(qXOpt, qYOpt, qTOpt)
      val queryXTY: Seq[Option[Extent[_]]] = Seq(qXOpt, qTOpt, qYOpt)

      def getRanges(name: String, curve: SpaceFillingCurve, query: Seq[Option[Extent[_]]]): List[IndexRange] = {
        val ranges: List[IndexRange] = curve.queryRanges(query).toList
        println(s"$name -> " + RichQueryRanges(ranges).toString)
        ranges
      }

      val ranges_Rtxy: List[IndexRange] = getRanges("Rtxy", Rtxy, queryTXY)
      ranges_Rtxy.size should be(22)

      val ranges_Rt_xy: List[IndexRange] = getRanges("Rt_xy", Rt_xy, queryTXY)
      ranges_Rt_xy.size should be(22)

      val ranges_Rxt_y: List[IndexRange] = getRanges("Rxt_y", Rxt_y, queryXTY)
      ranges_Rxt_y.size should be(22)

      val ranges_Rxy_t: List[IndexRange] = getRanges("Rxy_t", Rxy_t, queryXYT)
      ranges_Rxy_t.size should be(2)

      1 should be(1)
    }
  }
}

object ImplicitCasts {
  implicit class RichGMIR(gmir: GapMergedIndexRange) {
    override def toString: String = if (gmir.lower == gmir.upper) gmir.lower.toString
      else s"${gmir.lower.toString}-${gmir.upper.toString}(${gmir.size})"
  }

  implicit class RichQueryRanges(queryRanges: List[IndexRange]) {
    require(queryRanges != null)

    val numCells: Long = queryRanges.map(_.size).sum

    override def toString: String = {
      s"${queryRanges.size} entries, $numCells cells, ${queryRanges.map(qr => RichGMIR(qr.asInstanceOf[GapMergedIndexRange]).toString).mkString("[", ", ", "]")}"
    }
  }
}