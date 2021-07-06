package org.locationtech.sfcurve

import org.locationtech.sfcurve.Dimensions._
import org.locationtech.sfcurve.Utilities.CartesianProductIterable
import org.scalatest.MustMatchers.convertToAnyMustWrapper
import org.scalatest.{FunSpec, Matchers}

import java.text.SimpleDateFormat
import java.util.Date
import scala.util.{Failure, Success, Try}

class ComposedCurveSpec extends FunSpec with Matchers {
  val dtf = new SimpleDateFormat("yyyy-MM-dd HH:mm")
  val dt0: Date = dtf.parse("1988-02-12 11:45")
  val dt1: Date = dtf.parse("1998-08-28 22:15")
  val dt2: Date = dtf.parse("2001-04-23 12:30")
  val dt3: Date = dtf.parse("2019-02-12 08:11")

  val BitsPerDimension: Int = 20
  val X: Longitude = Longitude(BitsPerDimension)
  val Y: Latitude = Latitude(BitsPerDimension)
  val T: ExampleEra = ExampleEra(BitsPerDimension)

  val qTOpt: Option[Extent[_]] = Option(Extent(dt0, dt3))
  val qXOpt: Option[Extent[_]] = Option(Extent(-79.0, -78.0))
  val qYOpt: Option[Extent[_]] = Option(Extent(38.0, 39.0))

  def R(children: Discretizor*): RowMajorCurve = RowMajorCurve(children:_*)

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
      println(ranges.mkString("[", "; ", "]"))
    }
  }
}

/**
  * This is a naive implementation of a row-major curve for an arbitrary number of
  * dimensions.  It's really only used for unit testing.
  *
  * @param subordinates the child discretizors (either Dimension or SpaceFillingCurve)
  *                     from which to build this row-major curve
  */
case class RowMajorCurve(subordinates: Discretizor*) extends SpaceFillingCurve with InMemoryRangeConsolidator {
  val children: Vector[Discretizor] = subordinates.toVector

  val cardinalities: Vector[Long] = children.map(_.cardinality).foldRight(Vector(1L))((childCard, acc) => {
    acc :+ (acc.last * childCard)
  }).reverse.drop(1)

  def fold(subordinates: Seq[Long]): Long = {
    subordinates.zip(cardinalities).foldLeft(0L)((acc, pair) => pair match {
      case (childCoord, childCardinality) => acc + childCoord * childCardinality
    })
  }

  def unfold(index: Long): Vector[Long] = {
    cardinalities.foldLeft((Vector[Long](), index))((accPair, childCard) => accPair match {
      case (acc, indexLeft) => (
        acc :+ (indexLeft / childCard),
        indexLeft % childCard
      )
    })._1
  }

  // this is a BAD implementation because it iterates over every possible cell in the query range
  def indexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long], hints: Option[RangeComputeHints] = None): Seq[IndexRange] = {
    val maxGap: Int = hints.map(_.get(GapMergedIndexRange.HintsKeyMapGap)).getOrElse(0).toString.toInt

    val querySeqs: Seq[Seq[_]] = lowerCorner.zip(upperCorner).map {
      case (min, max) => (min to max).toSeq
    }

    val rawRanges: Iterator[IndexRange] = CartesianProductIterable(querySeqs).iterator.map(seq => {
      val coord: Seq[Long] = seq.map(_.asInstanceOf[Long])
      val index = fold(coord)
      GapMergedIndexRange(index, index, maxGap)
    })
    consolidateRanges(rawRanges, maxGap).toList
  }
}
