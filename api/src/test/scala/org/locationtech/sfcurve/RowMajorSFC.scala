package org.locationtech.sfcurve

import org.locationtech.sfcurve.Dimensions.{Discretizor, SpaceFillingCurve}
import org.locationtech.sfcurve.Utilities.CartesianProductIterable

/**
  * This is a naive implementation of a row-major curve for an arbitrary number of
  * dimensions.  It's only used for unit testing.
  *
  * Although it's a reasonable example of how to create a new SFC from the base
  * class, you should NOT use this implementation.  Instead, come up with a more
  * efficient range-query planner first, and swap that in.
  *
  * @param children the child discretizors (either Dimension or SpaceFillingCurve)
  *                 from which to build this row-major curve
  * @param maxGap the maximum allowable gap between query ranges
  */
case class RowMajorSFC(children: Vector[Discretizor], maxGap: Long = 0) extends SpaceFillingCurve with InMemoryRangeConsolidator {
  val ascendingCardinalities: Seq[Long] = children.map(_.cardinality).drop(1).scan(1L)(_ * _).reverse
  require(ascendingCardinalities.size == children.size)

  def fold(subordinates: Seq[Long]): Long = {
    // extra checks, as we're in unit tests
    subordinates.zip(children).foreach {
      case (index, discretizor) =>
        require(index >= 0)
        require(index < discretizor.cardinality)
    }
    // do the arithmetic
    subordinates.zip(ascendingCardinalities).map {
      case (idx, factor) => idx * factor
    }.sum
  }

  def unfold(index: Long): Vector[Long] = {
    ascendingCardinalities.foldLeft((index, Vector[Long]()))((t, ascCard) => t match {
      case (indexRemaining, subsSoFar) =>
        val idx: Long = indexRemaining / ascCard
        val rem: Long = indexRemaining % ascCard
        (rem, subsSoFar :+ idx)
    })._2
  }

  // this is a BAD implementation, because it's just for testing and needn't be efficient
  def indexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long], hints: Option[RangeComputeHints] = None): Seq[IndexRange] = {
    require(lowerCorner.size == arity)
    require(upperCorner.size == arity)

    val seqs: Seq[(Long, Long)] = lowerCorner.zip(upperCorner)
    val rangesIn: Seq[Seq[Long]] = seqs.map {
      case (a, b) => a to b
    }
    val coords: Iterator[Seq[Long]] = CartesianProductIterable(rangesIn).iterator.map(c => c.asInstanceOf[Seq[Long]])
    val indexes: Seq[IndexRange] = coords.map(coord => fold(coord)).toSeq.map(i => CoveredRange(i, i))
    val indexRanges = consolidateRanges(indexes.iterator, maxGap)
    indexRanges.toSeq
  }
}
