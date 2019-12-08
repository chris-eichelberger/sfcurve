package org.locationtech.sfcurve

import org.locationtech.sfcurve.Dimensions.SpaceFillingCurve
import org.locationtech.sfcurve.IndexRangeTree._

/**
  * This is the class responsible for supporting the in-memory
  * consolidation of index ranges into contiguous subsets.
  */
case class IndexRangeTree(indexRangeOpt: Option[GapMergedIndexRange] = None, leftChildOpt: Option[IndexRangeTree] = None, rightChildOpt: Option[IndexRangeTree] = None) {
  def getFarMin: Option[Long] = {
    leftChildOpt.map(_.getFarMin).getOrElse(indexRangeOpt.map(_.lower))
  }

  def getFarMax: Option[Long] = {
    rightChildOpt.map(_.getFarMax).getOrElse(indexRangeOpt.map(_.upper))
  }

  def add(gapMergedIndexRange: GapMergedIndexRange): IndexRangeTree = {
    if (indexRangeOpt.isEmpty) {
      return IndexRangeTree(gapMergedIndexRange, None, None)
    }

    val currentIndexRange = indexRangeOpt.get

    // already contained here
    if (currentIndexRange.contains(gapMergedIndexRange)) {
      return this
    }

    // inverse containment
    if (gapMergedIndexRange.contains(currentIndexRange)) {
      return IndexRangeTree(gapMergedIndexRange, leftChildOpt, rightChildOpt)
    }

    // overlaps
    if (currentIndexRange.overlaps(gapMergedIndexRange)) {
      val merged = currentIndexRange.merge(gapMergedIndexRange)
      return IndexRangeTree(merged, leftChildOpt, rightChildOpt)
    }

    // left
    if (gapMergedIndexRange.lt(currentIndexRange)) {
      return this.copy(leftChildOpt = leftChildOpt.map(_.add(gapMergedIndexRange)).orElse(IndexRangeTree(gapMergedIndexRange)))
    }

    // right
    if (gapMergedIndexRange.gt(currentIndexRange)) {
      return this.copy(rightChildOpt = rightChildOpt.map(_.add(gapMergedIndexRange)).orElse(IndexRangeTree(gapMergedIndexRange)))
    }

    throw new Exception(s"Should never reach here")
  }

  def add(tree: IndexRangeTree): IndexRangeTree = {
    tree.inOrderExtentIterator.foldLeft(this)((acc, extent) => {
      acc.add(extent)
    })
  }

  def inOrderExtentIterator: Iterator[GapMergedIndexRange] = {
    leftChildOpt.map(_.inOrderExtentIterator).getOrElse(EmptyExtentIterator) ++
      indexRangeOpt.map(_.iterator).getOrElse(EmptyExtentIterator) ++
      rightChildOpt.map(_.inOrderExtentIterator).getOrElse(EmptyExtentIterator)
  }

  def balance: IndexRangeTree = {
    val items = inOrderExtentIterator.toVector
    val indexes = getBalancedIndexOrdering(0, items.size - 1)

    indexes.foldLeft(IndexRangeTree())((acc, index) => {
      acc.add(items(index))
    })
  }
}

object IndexRangeTree {
  val EmptyExtentIterator: Iterator[GapMergedIndexRange] = Seq[GapMergedIndexRange]().iterator

  implicit def treeOption(tree: IndexRangeTree): Option[IndexRangeTree] = Option(tree)

  implicit def rangeOption(range: GapMergedIndexRange): Option[GapMergedIndexRange] = Option(range)

  def getBalancedIndexOrdering(lower: Int, upper: Int): Seq[Int] = {
    require(upper >= lower)
    if (upper == lower) return Seq(upper)
    val mid = (upper + lower) >> 1

    Nil ++
      Seq(mid) ++
      (if (mid > lower) getBalancedIndexOrdering(lower, mid - 1) else Nil) ++
      (if (mid < upper) getBalancedIndexOrdering(mid + 1, upper) else Nil)
  }
}


object RangeConsolidator {
  val EmptyRangeIndexIterator: Iterator[IndexRange] = List[IndexRange]().iterator
}

trait RangeConsolidator {
  // this trait only applies within a space-filling curve
  this: SpaceFillingCurve =>

  // given a list of ranges from each child, map them into ranges in your
  // own index space, and then post-process those ranges to ensure that we
  // have removed all gaps <= maxGap
  def consolidateRanges(subRanges: Iterator[IndexRange], maxGap: Long = 0): Iterator[IndexRange]
}

/**
  * Most of our base space-filling curves will have a query-range-planner that
  * knows how to emit ranges in order.  For that case, we can offer a simple
  * consolidator that does nothing but to return the identity.
  */
trait IdentityRangeConsolidator extends RangeConsolidator {
  // this trait only applies within a space-filling curve
  this: SpaceFillingCurve =>

  def consolidateRanges(subRanges: Iterator[IndexRange], maxGap: Long = 0): Iterator[IndexRange] = subRanges
}

trait InMemoryRangeConsolidator extends RangeConsolidator {
  // this trait only applies within a space-filling curve
  this: SpaceFillingCurve =>

  /**
    * To have gotten here, we have already retrieved all of the consolidated
    * index ranges from all children.  For example, given these raw values...
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
    * to the number of cells in the consolidated, output ranges.
    * That is:  (x:[1-5, 10-19], y:[100-119]) -> [z1-z2], [z3-z4], ..., [zN-zM]
    * such that sum of all the z-ranges is |x||y| = 15x20 = 300.
    * 2.  The maximum index from any child cannot exceed the maximum
    * possible index in the parent.  (The indexes *reported* by
    * children can exceed those *reported* by the parent, but only
    * when the child indexes remain smaller than the maximum
    * cardinality of the parent.)
    *
    * NB:  This default implementation handles everything IN MEMORY, which is
    * likely a _terrible_ way to do things (for SFCs of non-trivial size).
    * Use this routine at your own risk!
    *
    * @param subRanges per-child list of index ranges; the Cartesian product
    *                  of all these lists is the raw set of contiguous regions
    *                  for which indexes are sought
    * @param maxGap the maximum allowable gap between this curve's indexe
    *               ranges that are returned; any gap of less than or equal
    *               to this many grid-cells will be merged
    * @return
    */
  def consolidateRanges(subRanges: Iterator[IndexRange], maxGap: Long = 0): Iterator[IndexRange] = {
    // in-memory structure to maintain the collection of all valid ranges
    var indexRangeTree = IndexRangeTree()

    // accumulate these index ranges, growing the existing ranges appropriately
    subRanges.foreach(range => indexRangeTree = indexRangeTree.add(GapMergedIndexRange(range.lower, range.upper, maxGap)))

    indexRangeTree.inOrderExtentIterator
  }
}
