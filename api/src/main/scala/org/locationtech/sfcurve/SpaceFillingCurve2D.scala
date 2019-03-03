/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve

import org.locationtech.sfcurve.Dimensions.{Cell, Dimension, DimensionLike, Discretizor, Extent, Latitude, Longitude, NonNegativeLatitude, NonNegativeLongitude, SpaceFillingCurve}

class RangeComputeHints extends java.util.HashMap[String, AnyRef]

sealed trait IndexRange {
  def lower: Long
  def upper: Long
  def contained: Boolean
  def tuple = (lower, upper, contained)
  def size: Long = upper - lower + 1
}

case class CoveredRange(lower: Long, upper: Long) extends IndexRange {
  val contained = true
}

case class OverlappingRange(lower: Long, upper: Long) extends IndexRange {
  val contained = false
}

case class GapMergedIndexRange(lower: Long, upper: Long, maxGap: Long) extends IndexRange {
  // this is a degenerate method
  // TODO figure out how to thread this value through, if needed
  def contained: Boolean = false

  override def toString(): String = s"GapMergedIndexRange([$lower, $upper], gap $maxGap)"

  def contains(that: GapMergedIndexRange): Boolean =
    (lower + maxGap) <= that.lower && (upper - maxGap) >= that.upper

  // NB:  "overlaps" here means really "overlaps OR ABUTS"
  def overlaps(that: GapMergedIndexRange): Boolean =
    (lower - maxGap - 1) <= that.upper && (upper + maxGap + 1) >= that.lower

  def disjoint(that: GapMergedIndexRange): Boolean = !overlaps(that)

  def lt(that: GapMergedIndexRange): Boolean =
    upper < (that.lower - maxGap)

  def gt(that: GapMergedIndexRange): Boolean =
    lower > (that.upper + maxGap)

  def merge(that: GapMergedIndexRange): GapMergedIndexRange = {
    if (disjoint(that)) println(s"WARN:  Merging not-overlapping GapMergedIndexRanges $this and $that at upper-gap $maxGap")
    GapMergedIndexRange(Math.min(lower, that.lower), Math.max(upper, that.upper), maxGap)
  }
}

object GapMergedIndexRange {
  val HintsKeyMapGap: String = "MaximumAllowableGapBetweenQueryRanges"
}

object IndexRange {
  trait IndexRangeOrdering extends Ordering[IndexRange] {
    override def compare(x: IndexRange, y: IndexRange): Int = {
      val c1 = x.lower.compareTo(y.lower)
      if(c1 != 0) return c1
      val c2 = x.upper.compareTo(y.upper)
      if(c2 != 0) return c2
      0
    }
  }

  implicit object IndexRangeIsOrdered extends IndexRangeOrdering

  def apply(l: Long, u: Long, contained: Boolean): IndexRange =
    if(contained) CoveredRange(l, u)
    else          OverlappingRange(l, u)
}

/**
  *
  * @param bitsPrecision the number of bits to use PER DIMENSION to determine cardinality;
  *                      curve_cardinality = num_dimensions * (1 << bitsPrecision)
  */
abstract class SpaceFillingCurve2D(bitsPrecision: Int) extends SpaceFillingCurve {
  // we assume that the cardinality is specified in terms of the number of bits precision
  val xDimension = Longitude(1L << bitsPrecision)
  val yDimension = Latitude(1L << bitsPrecision)
  val children: Vector[Discretizor] = Vector(xDimension, yDimension)

  @deprecated("use 'index' instead", "SFCurve 2.0")
  def toIndex(x: Double, y: Double): Long = fold(Seq(xDimension.toBin(x), yDimension.toBin(y)))

  @deprecated("use 'inverseIndex' instead", "SFCurve 2.0")
  def toPoint(i: Long): (Double, Double) = {
    val cell = inverseIndex(i)
    (
      cell.extents.headOption.map(a => 0.5 * (a.min.asInstanceOf[Double] + a.max.asInstanceOf[Double])).get,
      cell.extents.lastOption.map(a => 0.5 * (a.min.asInstanceOf[Double] + a.max.asInstanceOf[Double])).get
    )
  }

  @deprecated("use 'indexRanges' instead", "SFCurve 2.0")
  def toRanges(xmin: Double, ymin: Double, xmax: Double, ymax: Double, hints: Option[RangeComputeHints] = None): Seq[IndexRange] = {
    val extXOpt = Option(Extent[Double](xmin, xmax))
    val extYOpt = Option(Extent[Double](ymin, ymax))
    queryRanges(Seq(extXOpt, extYOpt), hints).toSeq
  }
}
