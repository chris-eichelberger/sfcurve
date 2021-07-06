package org.locationtech.sfcurve

import org.locationtech.sfcurve.Dimensions._
import org.locationtech.sfcurve.Utilities.CartesianProductIterable
import org.scalatest.MustMatchers.convertToAnyMustWrapper
import org.scalatest.{FunSpec, Matchers}

import java.text.SimpleDateFormat
import java.util.Date
import scala.util.{Failure, Success, Try}

import ImplicitCasts._

class ComposedCurveSpec extends FunSpec with Matchers {
  val dtf = new SimpleDateFormat("yyyy-MM-dd HH:mm")
  val dt0: Date = dtf.parse("1988-02-12 11:45")
  val dt1: Date = dtf.parse("1998-08-28 22:15")
  val dt2: Date = dtf.parse("2001-04-23 12:30")
  val dt3: Date = dtf.parse("2019-02-12 08:11")

  val BitsPerDimension: Int = 10
  val CardinalityPerDimension: Long = 1L << BitsPerDimension.toLong
  val X: Longitude = Longitude(CardinalityPerDimension)
  val Y: Latitude = Latitude(CardinalityPerDimension)
  val T: ExampleEra = ExampleEra(CardinalityPerDimension)

  val qTOpt: Option[Extent[_]] = Option(Extent(dt0, dt3))
  val qXOpt: Option[Extent[_]] = Option(Extent(-79.0, -78.0))
  val qYOpt: Option[Extent[_]] = Option(Extent(38.0, 39.0))

  def R(children: Discretizor*): RowMajorSFC = RowMajorSFC(Vector(children:_*), 0)

  describe("???") {
    val Rtxy: SpaceFillingCurve = R(T, X, Y)
    val RZxy_t: SpaceFillingCurve = R(R(X, Y), T)
    val RtZxy: SpaceFillingCurve = R(T, R(X, Y))

    it("must have been able to instantiate the curves") {
      RtZxy must not be(null)
      RZxy_t must not be(null)
      Rtxy must not be(null)
    }

    it("must produce valid query ranges for different curves") {
      val ranges: List[IndexRange] = Rtxy.queryRanges(Seq(qTOpt, qXOpt, qYOpt)).toList
      println(RichQueryRanges(ranges).toString)
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