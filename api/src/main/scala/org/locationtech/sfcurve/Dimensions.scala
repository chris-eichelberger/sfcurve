package org.locationtech.sfcurve

import java.util.Date

object Dimensions {
  case class Extent[T](min: T, max: T, incMin: Boolean = true, incMax: Boolean = false)(implicit ev: Ordering[T]) {
    def contains(v: Any): Boolean = {
      val value: T = v.asInstanceOf[T]
      if (ev.lt(value, min)) return false
      if (ev.gt(value, max)) return false
      if (ev.eq(value, min) && !incMin) return false
      if (ev.eq(value, max) && !incMax) return false
      true
    }

    // ignores openness
    def isPoint: Boolean = ev.equiv(min, max)
  }

  trait DimensionLike[T] {
    def toBin(value: T, extent: Extent[T], cardinality: Long): Long
    def toExtent(idx: Long, extent: Extent[T], cardinality: Long): Extent[T]
  }

  implicit object LongDimensionLike extends DimensionLike[Long] {
    def toBin(value: Long, extent: Extent[Long], cardinality: Long): Long = {
      require(extent.contains(value))
      Math.round((value - extent.min) * cardinality.toDouble / (extent.max - extent.min))
    }

    def toExtent(bin: Long, extent: Extent[Long], cardinality: Long): Extent[Long] = {
      require(bin >= 0 && bin < cardinality)
      val x0: Long = extent.min + Math.floor(bin.toDouble / cardinality.toDouble * (extent.max - extent.min)).toLong
      val x1: Long = extent.min + Math.floor((1.0 + bin.toDouble) / cardinality.toDouble * (extent.max - extent.min)).toLong
      Extent(x0, x1, incMin = true, incMax = false)
    }
  }

  implicit object DoubleDimensionLike extends DimensionLike[Double] {
    def toBin(value: Double, extent: Extent[Double], cardinality: Long): Long = {
      require(extent.contains(value))
      Math.round((value - extent.min) * cardinality.toDouble / (extent.max - extent.min))
    }

    def toExtent(bin: Long, extent: Extent[Double], cardinality: Long): Extent[Double] = {
      require(bin >= 0 && bin < cardinality)
      val x0: Double = extent.min + bin.toDouble / cardinality.toDouble * (extent.max - extent.min)
      val x1: Double = extent.min + (1.0 + bin.toDouble) / cardinality.toDouble * (extent.max - extent.min)
      Extent(x0, x1, incMin = true, incMax = false)
    }
  }

  implicit def date2long(date: Date): Long = date.getTime
  implicit def long2date(long: Long): Date = new Date(long)
  implicit def dex2lex(extent: Extent[Date]): Extent[Long] =
    Extent[Long](date2long(extent.min), date2long(extent.max), extent.incMin, extent.incMax)
  implicit def lex2dex(extent: Extent[Long]): Extent[Date] =
    Extent[Date](long2date(extent.min), long2date(extent.max), extent.incMin, extent.incMax)

  implicit object DateDimensionLike extends DimensionLike[Date] {
    def toBin(value: Date, extent: Extent[Date], cardinality: Long): Long =
      LongDimensionLike.toBin(date2long(value), extent, cardinality)

    def toExtent(bin: Long, extent: Extent[Date], cardinality: Long): Extent[Date] =
      LongDimensionLike.toExtent(bin, extent, cardinality)
  }

  trait Discretizor {
    def cardinality: Long
    def arity: Int
    def index(values: Seq[DimensionLike[_]]): Long
    def inverseIndex(index: Long): Cell
    def getQueryRanges(extents: Seq[Option[Extent[_]]], hints: Option[RangeComputeHints] = None): Seq[IndexRange]
  }

  /**
    * This is a discretizor that knows how to partition a typed
    * space into bins (and backwards).  Examples:  Double-valued
    * spaces such as real-number ranges; Date-valued spaces; etc.
    *
    * @tparam T describes the space being partitioned
    */
  trait Dimension[T] extends Discretizor {
    def ev: DimensionLike[T]
    def extent: Extent[T]
    def fullIndexRange: IndexRange = CoveredRange(0, cardinality - 1)

    val arity: Int = 1

    def toBin(value: T): Long = ev.toBin(value, extent, cardinality)

    def toExtent(bin: Long): Extent[T] = ev.toExtent(bin, extent, cardinality)

    def index(values: Seq[DimensionLike[_]]): Long = {
      require(values.size == 1)
      values.head match {
        case value: T => toBin(value)
        case _ => throw new Exception("Invalid value type for index")
      }
    }

    def inverseIndex(index: Long): Cell = {
      val extent: Extent[Dimension[_]] = toExtent(index).asInstanceOf[Extent[Dimension[_]]]
      val extents: Vector[Extent[Dimension[_]]] = Vector(extent)
      Cell(extents)
    }

    // dummy identity function; some dimensions may override this
    def normalize(value: T): T = value

    // a dimension will always return a contiguous range of bin indexes from a
    // contiguous range of raw input values
    def getQueryRanges(extents: Seq[Option[Extent[_]]], hints: Option[RangeComputeHints] = None): Seq[IndexRange] = {
      require(extents.size == 1)
      val subExtent = extents.head
      if (subExtent.isEmpty) return Seq(fullIndexRange)
      if (subExtent.contains(extent)) return Seq(fullIndexRange)
      Seq(CoveredRange(toBin(subExtent.get.min.asInstanceOf[T]), toBin(subExtent.get.max.asInstanceOf[T])))
    }
  }

  // defined over [-180.0, 180.0)
  case class Longitude(cardinality: Long) extends Dimension[Double] {
    val ev: DimensionLike[Double] = implicitly[DimensionLike[Double]]
    val extent: Extent[Double] = Extent[Double](-180.0, 180.0, incMin = true, incMax = false)
    override def normalize(value: Double): Double = {
      var v = value
      while (v < -180.0) v = 360.0
      while (v >= 180.0) v -= 360.0
      v
    }
  }

  // defined over [-90.0, 90.0)
  case class Latitude(cardinality: Long) extends Dimension[Double] {
    val ev: DimensionLike[Double] = implicitly[DimensionLike[Double]]
    val extent: Extent[Double] = Extent[Double](-180.0, 180.0, incMin = true, incMax = false)
    override def normalize(value: Double): Double = {
      var v = value
      while (v < -90.0) v = 180.0
      while (v >= 90.0) v -= 180.0
      v
    }
  }

  case class Cell(extents: Vector[Extent[Dimension[_]]] = Vector()) {
    val numDimensions: Int = extents.size

    def contains(point: Vector[Dimension[_]]): Boolean = {
      if (point.size != numDimensions) return false
      extents.zip(point).forall {
        case (extent, coord) => extent.contains(coord)
      }
    }

    def +(that: Cell): Cell = Cell(extents ++ that.extents)
  }

  trait SpaceFillingCurve extends Discretizor {
    def children: Vector[Discretizor]
    require(children != null && children.nonEmpty)

    // how many possible values can the `index` method return?
    def cardinality: Long = children.map(_.cardinality).product

    // how many dimension leaf-nodes are there?
    def arity: Int = children.map(_.arity).sum

    // these routines are the heart of the SFC;
    // they will be made concrete per descendant curve
    def fold(subordinates: Seq[Long]): Long
    def unfold(index: Long): Vector[Long]
    def getIndexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long]): Seq[IndexRange]

    def index(values: Seq[DimensionLike[_]]): Long = {
      require(values.size == arity)

      // defer to your children, be they SFCs or Dimensions
      val subordinates: Seq[Long] = children.foldLeft((values, Seq[Long]()))((acc, child) => acc match {
        case (valuesLeft: Seq[DimensionLike[_]], subsSoFar: Seq[Long]) =>
          val childValues: Seq[DimensionLike[_]] = valuesLeft.take(child.arity)
          val childIndex: Long = child.index(childValues)
          (valuesLeft.drop(child.arity), subsSoFar :+ childIndex)
      })._2

      // roll up the child indexes into a single index
      fold(subordinates)
    }

    def inverseIndex(index: Long): Cell = {
      // unfold this index to get the child indexes to recurse
      val subordinates: Seq[Long] = unfold(index)

      // this should always be true
      require(subordinates.size == children.size)

      val cell: Cell = children.zip(subordinates).foldLeft(Cell())((acc, t) => t match {
        case (child, childIndex) =>
          acc + child.inverseIndex(childIndex)
      })

      // sanity check
      require(cell.extents.size == arity)

      cell
    }

    def getQueryRanges(extents: Seq[Option[Extent[_]]], hints: Option[RangeComputeHints] = None): Seq[IndexRange] = {
      // assume that unqueried dimensions are either:
      // 1.  an option set to its full extent; or
      // 2.  a None
      require(extents.size == arity)

      // defer to your children, be they SFCs or Dimensions
      val subordinates: Seq[Seq[IndexRange]] = children.foldLeft((extents, Seq[Seq[IndexRange]]()))((acc, child) => acc match {
        case (extentsLeft: Seq[Option[Extent[_]]], subsSoFar: Seq[IndexRange]) =>
          val childExtents: Seq[Option[Extent[_]]] = extentsLeft.take(child.arity)
          val childRanges: Seq[IndexRange] = child.getQueryRanges(childExtents)
          (extentsLeft.drop(child.arity), subsSoFar :+ childRanges)
      })._2

      // roll up the child index-ranges into a single index-range
      consolidateRanges(subordinates)
    }

    // copied from https://github.com/aheyne/geomesa/blob/master/geomesa-utils/src/main/scala/org/locationtech/geomesa/utils/iterators/CartesianProductIterable.scala
    /**
      * Can create an iterator over all combinations of items from a list-of-lists.
      * Because the final list of combinations can be large, we allow for a safe
      * way to query the list size that is independent of the iterator itself.
      * (That is, asking for the size does not exhaust any iterator.)
      *
      * NB:  The first sequence is the least significant; that is, it will
      * increment fast while the last sequence is the most significant (will
      * increment slowly).
      *
      * @param seqs the list-of-lists whose items are to be recombined
      * @tparam T the type of items
      */
    case class CartesianProductIterable(seqs: Seq[Seq[_]]) extends Iterable[Seq[_]] {
      lazy val expectedSize: Long = seqs.map(_.size.toLong).product

      def iterator: Iterator[Seq[_]] = new Iterator[Seq[_]] {
        val n: Int = seqs.size
        val maxes: Vector[Int] = seqs.map(seq => seq.size).toVector
        val indexes = new scala.collection.mutable.ArraySeq[Int](seqs.size)
        var nextItem: Seq[_] = if (isValid) realize else null

        def isValid: Boolean = (0 until n).forall(i => indexes(i) < maxes(i))

        def realize: Seq[_] = (0 until n).map(i => seqs(i)(indexes(i)))

        def hasNext: Boolean = nextItem != null

        def next(): Seq[_] = {
          if (nextItem == null) throw new Exception("Iterator exhausted")
          val result = nextItem

          // advance the internal state
          nextItem = null
          var j = 0
          var done = false
          while (j < n && !done) {
            indexes(j) = indexes(j) + 1
            if (indexes(j) >= maxes(j)) {
              indexes(j) = 0
              j = j + 1
            } else {
              done = true
            }
          }
          if (done || j < n) nextItem = realize

          result
        }
      }
    }

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
      * 1.  The total number of cells in the input ranges MUST be equal
      *     to the number of cells in the consolidated, output ranges.
      * 2.  The maximum index from any child cannot exceed the maximum
      *     possible index in the parent.  (The indexes *reported* by
      *     children can exceed those *reported* by the parent, but only
      *     when the child indexes remain smaller than the maximum
      *     cardinality of the parent.)
      *
      * @param subRanges
      * @return
      */
    def consolidateRanges(subRanges: Seq[Seq[IndexRange]]): Seq[IndexRange] = {
      // there must be one sequence, even if empty, per child
      require(subRanges.size == children.size)

      val productItr: Iterator[Seq[_]] = CartesianProductIterable(subRanges).iterator
      for (combination <- productItr) {
        // find the minimum index
        val lowerCorner = combination.map(dimRange => dimRange.asInstanceOf[IndexRange].lower)
        val upperCorner = combination.map(dimRange => dimRange.asInstanceOf[IndexRange].upper)
        val ranges: Seq[IndexRange] = getIndexRanges(lowerCorner, upperCorner)
      }
    }
  }
}
