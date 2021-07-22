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

/**
  * Represents a two-dimensional Hilbert curve.  This class is just a facade, with the
  * Uzaygezen library doing all of the heavy lifting.
  *
  * NB:  Because of the parent contract, this will always be a square Hilbert curve.
  *
  * @param bitsPerDimension the number of bits resolution to use in EACH of the two dimensions
  */
class HilbertCurve2D(bitsPerDimension: Int) extends SpaceFillingCurve2D(bitsPerDimension) with IdentityRangeConsolidator {
  require(bitsPerDimension > 0, "Bits/dimension is negative; this probably means you tried to pass in a Long value that was too large")
  require(bitsPerDimension <= 31, s"Bits/dimension is too large; must be no more than 31")
  val chc = new CompactHilbertCurve(Array(bitsPerDimension, bitsPerDimension))

  override val name: String = "Hilbert"

  def fold(subordinates: Seq[Long]): Long = {
    require(subordinates.length == 2)

    val normX = subordinates.head
    val normY = subordinates.last

    val p =
      Array[BitVector](
        BitVectorFactories.OPTIMAL(bitsPerDimension),
        BitVectorFactories.OPTIMAL(bitsPerDimension)
      )

    p(0).copyFrom(normX)
    p(1).copyFrom(normY)

    val hilbert = BitVectorFactories.OPTIMAL.apply(bitsPerDimension * 2)

    chc.index(p,0,hilbert)
    hilbert.toLong
  }

  def unfold(index: Long): Vector[Long] = {
    val h = BitVectorFactories.OPTIMAL.apply(bitsPerDimension*2)
    h.copyFrom(index)
    val p =
      Array[BitVector](
        BitVectorFactories.OPTIMAL(bitsPerDimension),
        BitVectorFactories.OPTIMAL(bitsPerDimension)
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

    val chc = new CompactHilbertCurve(Array[Int](bitsPerDimension, bitsPerDimension))
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
