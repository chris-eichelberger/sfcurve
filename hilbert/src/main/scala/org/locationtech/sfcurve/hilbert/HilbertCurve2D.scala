/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve.hilbert

import java.io.Serializable

import com.google.common.base.{Function, Functions}
import com.google.common.collect.ImmutableList
import com.google.uzaygezen.core._
import com.google.uzaygezen.core.ranges._
import org.locationtech.sfcurve._

class HilbertCurve2DProvider extends SpaceFillingCurveProvider {
  override def canProvide(name: String): Boolean = "hilbert" == name

  override def build2DSFC(args: Map[String, Serializable]): SpaceFillingCurve2D =
    new HilbertCurve2D(args(HilbertCurve2DProvider.RESOLUTION_PARAM).asInstanceOf[Int])
}

object HilbertCurve2DProvider {
  val RESOLUTION_PARAM = "hilbert.resolution"
}

class HilbertCurve2D(resolution: Int) extends SpaceFillingCurve2D(Dimensions.bitsFromCardinality(resolution).toInt) with IdentityRangeConsolidator {
  require(resolution > 0, "Resolution is negative; this probably means you tried to pass in a Long value that was too large")
  require(resolution <= (Int.MaxValue >> 1), s"Resolution is too large; must be no more than ${Int.MaxValue >> 1}")
  val bitsPerDimension: Int = Dimensions.bitsFromCardinality(resolution).toInt
  val precision: Long = Dimensions.bitsFromCardinality(resolution)
  val chc = new CompactHilbertCurve(Array(bitsPerDimension, bitsPerDimension))

  val name: String = "Hilbert"

  def fold(subordinates: Seq[Long]): Long = {
    require(subordinates.length == 2)

    val normX = subordinates.head
    val normY = subordinates.last

    val p =
      Array[BitVector](
        BitVectorFactories.OPTIMAL(resolution),
        BitVectorFactories.OPTIMAL(resolution)
      )

    p(0).copyFrom(normX)
    p(1).copyFrom(normY)

    val hilbert = BitVectorFactories.OPTIMAL.apply(resolution * 2)

    chc.index(p,0,hilbert)
    hilbert.toLong
  }

  def unfold(index: Long): Vector[Long] = {
    val h = BitVectorFactories.OPTIMAL.apply(resolution*2)
    h.copyFrom(index)
    val p =
      Array[BitVector](
        BitVectorFactories.OPTIMAL(resolution),
        BitVectorFactories.OPTIMAL(resolution)
      )

    chc.indexInverse(h,p)

    Vector(p(0).toLong, p(1).toLong)
  }

  // the Uzaygezen implementation returns ranges in consistently /descending/ order
  def indexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long], hints: Option[RangeComputeHints] = None): Seq[IndexRange] = {
    require(lowerCorner.size == 2)
    val minNormalizedLongitude = lowerCorner.head
    val minNormalizedLatitude = lowerCorner.last
    require(upperCorner.size == 2)
    val maxNormalizedLongitude = upperCorner.head
    val maxNormalizedLatitude = upperCorner.last

    // dummy check for a single-cell overlap (which Uzaygezen does /not/ like)
    if (minNormalizedLongitude == maxNormalizedLongitude && minNormalizedLatitude == maxNormalizedLatitude) {
      val idx = fold(Seq(minNormalizedLongitude, minNormalizedLatitude))
      return Seq(IndexRange(idx, idx, contained = true))
    }

    // dummy check for a single-cell strip (which Uzaygezen does /not/ like)
    if (minNormalizedLongitude == maxNormalizedLongitude) {
      // vertical strip
      return (for (normalizedLat <- minNormalizedLatitude to maxNormalizedLatitude) yield {
        val idx = fold(Seq(minNormalizedLongitude, normalizedLat))
        IndexRange(idx, idx, contained = true)
      }).toSeq
    }
    if (minNormalizedLatitude == maxNormalizedLatitude) {
      // horizontal strip
      return (for (normalizedLon <- minNormalizedLongitude to maxNormalizedLongitude) yield {
        val idx = fold(Seq(normalizedLon, minNormalizedLatitude))
        IndexRange(idx, idx, contained = true)
      }).toSeq
    }
    if (minNormalizedLongitude == maxNormalizedLongitude || minNormalizedLatitude == maxNormalizedLatitude) {
      throw new Exception(s"Single-cell strip should not be possible:  long ($minNormalizedLongitude, $maxNormalizedLongitude), lat ($minNormalizedLatitude, $maxNormalizedLatitude)")
    }

    val chc = new CompactHilbertCurve(Array[Int](precision.toInt, precision.toInt))
    val region = new java.util.ArrayList[LongRange]()

    region.add(LongRange.of(minNormalizedLongitude,maxNormalizedLongitude))
    region.add(LongRange.of(minNormalizedLatitude,maxNormalizedLatitude))

    val zero = new LongContent(0L)
    val LongRangeIDFunction: Function[LongRange, LongRange] = Functions.identity()

    val inspector =
      SimpleRegionInspector.create(
        ImmutableList.of(region),
        new LongContent(1L),
        LongRangeIDFunction,
        LongRangeHome.INSTANCE,
        zero
      )

    val combiner =
      new PlainFilterCombiner[LongRange, java.lang.Long, LongContent, LongRange](LongRange.of(0, 1))

    val queryBuilder = BacktrackingQueryBuilder.create(inspector, combiner, Int.MaxValue, true, LongRangeHome.INSTANCE, zero)

    chc.accept(new ZoomingSpaceVisitorAdapter(chc, queryBuilder))

    val query = queryBuilder.get()

    val ranges = query.getFilteredIndexRanges

    //result
    var result = List[IndexRange]()
    val itr = ranges.iterator

    while(itr.hasNext) {
      val l = itr.next()
      val range = l.getIndexRange
      val start = range.getStart.asInstanceOf[Long]
      val end   = range.getEnd.asInstanceOf[Long]
      val contained = l.isPotentialOverSelectivity
      result = IndexRange(start, end, contained) :: result
    }
    result
  }

}
