package org.locationtech.sfcurve

import org.locationtech.sfcurve.Dimensions.{Discretizor, SpaceFillingCurve}
import org.locationtech.sfcurve.Utilities.CartesianProductIterable

import scala.collection.immutable.TreeMap

/**
  * This is the class responsible for supporting the in-memory
  * consolidation of index ranges into contiguous subsets.
  *
  * @param byFirst
  * @param byLast
  */
case class IndexedRanges(maxGap: Long, byFirst: TreeMap[Long, IndexRange], byLast: TreeMap[Long, IndexRange]) extends Iterable[IndexRange] {
  def +(indexRange: IndexRange): IndexedRanges = {
    
  }

  def iterator: Iterator[IndexRange] = byFirst.valuesIterator
}

trait RangeConsolidator {
  this: SpaceFillingCurve =>
  def children: Vector[Discretizor]
  def getIndexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long]): Seq[IndexRange]
  def consolidateRanges(subRanges: Seq[Seq[IndexRange]]): Seq[IndexRange]
}

trait InMemoryRangeConsolidator extends RangeConsolidator {
  /**
    * To have gotten here, we have already retrieved all of the consolidated
    * index ranges from all children.  For example, given these raw values:
    *
    *   Z(
    *     H(x: [-80, -79], y: [38, 39]),  // curve
    *     t: [2018-12-01, 2018-12-31]     // dimension
    *   )
    *
    * it might produce these index ranges to be consolidated by the parent
    * Z-order curve (values are made up for example purposes):
    *
    *   Z(
    *     H:  1001-1003, 1021-1023, 2017, 2018, 3001  // curves can return disjoint ranges
    *     t:  65530-65535                             // dimensions return only zero or one contiguous range
    *   )
    *
    * The *naive* process the Z-order curve has to do next is to consider
    * the full cartesian product of the start/stop pairs from its children,
    * and compute the ranges that fall inside, consolidating them as
    * appropriate.  For this case (values are made up for this example):
    *
    *   order_and_consolidate(
    *     Z(1001, 65530) .. Z(1003, 65535) -> 1048572-1048575, 1048604, 1048712-1048800
    *     Z(1021, 65530) .. Z(1023, 65535) -> 11-15, 19-31, 90
    *     Z(2017, 65530) .. Z(2018, 65535) -> 1048605-1048610, 1048620-1048630
    *     Z(3001, 65530) .. Z(3001, 65535) -> 1048619, 91-93
    *   )
    *
    * would consolidate to this final set of (made up, example!) ranges:
    *
    *   11-15
    *   19-31
    *   90-93
    *   1048572-1048575
    *   1048604-1048610
    *   1048619-1048630
    *   1048712-1048800
    *
    * This example, by virtue of being made up from whole cloth, does not honor
    * the constraints that you would expect in a real curve.  To wit:
    *
    * 1.  The product of the number of cells in the input ranges MUST be equal
    *     to the number of cells in the consolidated, output ranges.
    *     That is:  (x:[1-5, 10-19], y:[100-119]) -> [z1-z2], [z3-z4], ..., [zN-zM]
    *     such that sum of all the z-ranges is |x||y| = 15x20 = 300.
    * 2.  The maximum index from any child cannot exceed the maximum
    *     possible index in the parent.  (The indexes *reported* by
    *     children can exceed those *reported* by the parent, but only
    *     when the child indexes remain smaller than the maximum
    *     cardinality of the parent.)
    *
    * NB:  This default implementation handles everything IN MEMORY, which is
    * likely a _terrible_ way to do things (for SFCs of non-trivial size).
    * Use this routine at your own risk!
    *
    * @param subRanges
    * @return
    */
  def consolidateRanges(subRanges: Seq[Seq[IndexRange]]): Seq[IndexRange] = {
    // there must be one sequence, per child
    require(subRanges.size == children.size)

    val numCellsExpected = subRanges.map {
      case r: Seq[IndexRange] if r.isEmpty => 0
      case r: Seq[IndexRange]              => r.map(iRange => iRange.size).sum
    }.product

    // degenerate case:  at least one dimension has no coverage, so the
    // intersection must be empty
    if (numCellsExpected == 0L) return Seq()

    val productItr: Iterator[Seq[_]] = CartesianProductIterable(subRanges).iterator
    for (combination <- productItr) {
      // find the query rectangle in (input) index space
      val lowerCorner = combination.map(dimRange => dimRange.asInstanceOf[IndexRange].lower)
      val upperCorner = combination.map(dimRange => dimRange.asInstanceOf[IndexRange].upper)

      // find the query ranges in (output) index space that underlie the query rectangle
      val ranges: Seq[IndexRange] = getIndexRanges(lowerCorner, upperCorner)

      // accumulate these index ranges, growing the existing ranges appropriately
    }
  }
}
