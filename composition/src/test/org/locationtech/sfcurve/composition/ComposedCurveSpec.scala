package org.locationtech.sfcurve.composition

import org.locationtech.sfcurve.Dimensions._
import org.locationtech.sfcurve.composition.ImplicitCasts._
import org.locationtech.sfcurve.hilbert.HilbertCurve2D
import org.locationtech.sfcurve.zorder.{ZCurve2D, ZCurve3D}
import org.locationtech.sfcurve.{Dimensions, IndexRange, RowMajorSFC}
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
  val dt4: Date = dtf.parse("1988-02-22 11:45")

  /*
  Three-dimensional curves will come in two shapes:
  1.  C(a, a, a)
  2.  C(C(b, b), 2b) or C(2b, C(b, b))

  Where "a" and "b" are bit counts.  When a=4, for example, the first form will
  have a top-level curve of 12 bits resolution (cardinality of 4096).  When b=3,
  similarly, the second form will have a top-level curve of 12 bits resolution
  (cardinality of 4096).

  The implementations of Z-order and Hilbert assume square curves.

  In order for the top-level cardinality of the two curve forms to be the same,
  we have to make sure that 3a = 4b.  For example, a=8 and b=6 would work,
  which would give both curve forms a total resolution of 24 bits.
  (To say it out loud, "a" has to be divisible by 4, so then "b = 3a/4".)
   */
  
  // for curve form 1:  C(a, a, a)
  val BitsPerDimensionA: Int = 12  // must be divisible by 4
  val CardinalityPerDimensionA: Long = 1L << BitsPerDimensionA
  val XA: Longitude = Longitude(CardinalityPerDimensionA)
  val YA: Latitude = Latitude(CardinalityPerDimensionA)
  val TA: DateDimension = ExampleEra(CardinalityPerDimensionA)

  // for curve form 2:  C(C(b, b), 2b) or C(2b, C(b, b))
  val BitsPerDimensionB: Int = BitsPerDimensionA * 3 / 4
  val CardinalityPerDimensionB: Long = 1L << BitsPerDimensionB
  val LargeBitsPerDimensionB: Int = BitsPerDimensionB << 1
  val LargeCardinalityPerDimensionB: Long = 1L << LargeBitsPerDimensionB
  val XB: Longitude = Longitude(CardinalityPerDimensionB)
  val YB: Latitude = Latitude(CardinalityPerDimensionB)
  val TB: DateDimension = ExampleEra(CardinalityPerDimensionB)
  val LargeXB: Longitude = Longitude(LargeCardinalityPerDimensionB)
  val LargeYB: Latitude = Latitude(LargeCardinalityPerDimensionB)
  val LargeTB: DateDimension = ExampleEra(LargeCardinalityPerDimensionB)

  // 38.029925749354696, -78.47741019445668 is approximately Charlottesville City Hall
  val qTOpt: Option[Extent[_]] = Option(Extent(dt0, dt4))
  val dx: Double = 0.05
  val xMid: Double = -78.52690829423837
  val xMin: Double = xMid - dx
  val xMax: Double = xMid + dx
  val dy: Double = dx
  val yMid: Double = 38.029925749354696
  val yMin: Double = yMid - dy
  val yMax: Double = yMid + dy
  val qXOpt: Option[Extent[_]] = Option(Extent(xMin, xMax))
  val qYOpt: Option[Extent[_]] = Option(Extent(yMin, yMax))

  def R(children: Discretizor*): RowMajorSFC = RowMajorSFC(Vector(children:_*))

  def getRanges(name: String, curve: SpaceFillingCurve, query: Seq[Option[Extent[_]]], optNumEntries: Option[Long] = None, optNumCells: Option[Long] = None, printCellDetails: Boolean = false): List[IndexRange] = {
    val ranges: List[IndexRange] = curve.queryRanges(query).toList
    println(s"$name -> " + RichQueryRanges(ranges).toString)

    if (printCellDetails) {
      for (ir <- ranges.take(ImplicitCasts.MaxRanges)) {
        print(s"  Index range:  ${RichIR(ir).toString}")
        if (ir.size == 1) {
          val cell = curve.inverseIndex(ir.lower)
          println(s" -> $cell")
        } else {
          println(":")
          for (idx <- (ir.lower to ir.upper)) {
            val cell = curve.inverseIndex(idx)
            println(s"    $idx -> $cell")
          }
        }
      }
    }

    optNumEntries match {
      case Some(numEntries) => ranges.size must be(numEntries)
      case None =>
    }

    optNumCells match {
      case Some(numCells) => ranges.numCells must be(numCells)
      case None =>
    }

    ranges
  }

  describe("composed-curve variants using only R, x, y, and t") {
    val Rtxy: SpaceFillingCurve = R(TA, XA, YA)
    val Rxy_t: SpaceFillingCurve = R(R(XB, YB), LargeTB)
    val Rt_xy: SpaceFillingCurve = R(LargeTB, R(XB, YB))
    val Rxt_y: SpaceFillingCurve = R(R(XB, TB), LargeYB)

    it("must have been able to instantiate the curves") {
      Rt_xy must not be null
      Rxy_t must not be null
      Rtxy must not be null
      Rxt_y must not be null
    }

    it("must produce valid query ranges for different curves") {
      val queryTXY: Seq[Option[Extent[_]]] = Seq(qTOpt, qXOpt, qYOpt)
      val queryXYT: Seq[Option[Extent[_]]] = Seq(qXOpt, qYOpt, qTOpt)

      // form A uses bigger cells, so there should be fewer
      // (but all A-form curves should have the same number of cells)
      getRanges("Rtxy", Rtxy, queryTXY, Option(9), Option(27))

      // form B uses smaller cells, so there should be more
      // (but all B-form curves should have the same number of cells)
      getRanges("Rt_xy", Rt_xy, queryTXY, Some(73), Some(73))
      getRanges("Rxy_t", Rxy_t, queryXYT, Some(1), Some(73))
    }
  }

  // only supports 2- and 3-dimensional curves
  def Z(children: Discretizor*): SpaceFillingCurve = {
    val subordinates: Vector[Discretizor] = Vector(children:_*)
    require(subordinates.size >= 2 && subordinates.size <= 3)
    subordinates.size match {
      case 2 =>
        new ZCurve2D(subordinates.head.cardinality.toInt) {
          override val children: Vector[Discretizor] = subordinates
        }
      case 3 =>
        new ZCurve3D(subordinates.head.cardinality.toInt) {
          override val children: Vector[Discretizor] = subordinates
        }
    }
  }

  // only supports 2- and 3-dimensional curves
  def H(children: Discretizor*): SpaceFillingCurve = {
    val subordinates: Vector[Discretizor] = Vector(children:_*)
    require(subordinates.size == 2)
    new HilbertCurve2D(Dimensions.bitsFromCardinality(subordinates.head.cardinality).toInt) {
      override val children: Vector[Discretizor] = subordinates
    }
  }

  describe("composed curves using R, Z, H, x, y, and t") {
    val Ztxy: SpaceFillingCurve = Z(TA, XA, YA)
    val RtZxy: SpaceFillingCurve = R(LargeTB, Z(XB, YB))
    val RtHxy: SpaceFillingCurve = R(LargeTB, H(XB, YB))
    val HtZxy: SpaceFillingCurve = H(LargeTB, Z(XB, YB))
    val RZxy_t: SpaceFillingCurve = R(Z(XB, YB), LargeTB)
    val RHxy_t: SpaceFillingCurve = R(H(XB, YB), LargeTB)
    val HZxy_t: SpaceFillingCurve = H(Z(XB, YB), LargeTB)

    it("must be able to instantiate the curves") {
      RtZxy must not be null
      RtHxy must not be null
      Ztxy must not be null
      HtZxy must not be null
      RZxy_t must not be null
      RHxy_t must not be null
      HZxy_t must not be null
    }

    it("must produce valid query ranges for different curves") {
      // background
      getRanges(s"t:$LargeCardinalityPerDimensionB", R(LargeTB), Seq(qTOpt), None, None, printCellDetails = true)
      getRanges(s"Zxy:${CardinalityPerDimensionB*CardinalityPerDimensionB}", Z(XB, YB), Seq(qXOpt, qYOpt), None, None, printCellDetails = true)
      getRanges(s"Hxy:${CardinalityPerDimensionB*CardinalityPerDimensionB}", H(XB, YB), Seq(qXOpt, qYOpt), None, None, printCellDetails = true)
      getRanges(s"t$CardinalityPerDimensionB", R(TB), Seq(qTOpt), None, None, printCellDetails = true)
      getRanges(s"x$CardinalityPerDimensionB", R(XB), Seq(qXOpt), None, None, printCellDetails = true)
      getRanges(s"y$CardinalityPerDimensionB", R(YB), Seq(qYOpt), None, None, printCellDetails = true)

      // form A should have fewer cells
      getRanges("Ztxy", Ztxy, Seq(qTOpt, qXOpt, qYOpt), Some(13), Some(27))

      // all B-form curves should have similar numbers of cells
      getRanges("RtZxy", RtZxy, Seq(qTOpt, qXOpt, qYOpt), Some(73), Some(73))
      getRanges("RtHxy", RtHxy, Seq(qTOpt, qXOpt, qYOpt), Some(73), Some(73))
      getRanges("HtZxy", HtZxy, Seq(qTOpt, qXOpt, qYOpt), Some(73), Some(73))
      getRanges("RZxy_t", RZxy_t, Seq(qXOpt, qYOpt, qTOpt), Some(1), Some(73))
      getRanges("RHxy_t", RHxy_t, Seq(qXOpt, qYOpt, qTOpt), Some(1), Some(73))
      getRanges("HZxy_t", HZxy_t, Seq(qXOpt, qYOpt, qTOpt), Some(73), Some(73))
    }
  }
}

object ImplicitCasts {
  val MaxRanges: Int = 15

  implicit class RichIR(ir: IndexRange) {
    override def toString: String = if (ir.lower == ir.upper) ir.lower.toString
      else s"${ir.lower.toString}-${ir.upper.toString}(${ir.size})"
  }

  implicit class RichQueryRanges(queryRanges: List[IndexRange]) {
    require(queryRanges != null)

    val numCells: Long = queryRanges.map(_.size).sum

    override def toString: String = {
      if (queryRanges.size <= MaxRanges) {
        // exhaustive
        s"${queryRanges.size} entries, $numCells cells, ${queryRanges.map(qr => RichIR(qr).toString).mkString("[", ", ", "]")}"
      } else {
        // summary only
        s"${queryRanges.size} entries, $numCells cells, [${queryRanges.take(MaxRanges).map(qr => RichIR(qr).toString).mkString(", ")}, ... elided...]"
      }
    }
  }
}