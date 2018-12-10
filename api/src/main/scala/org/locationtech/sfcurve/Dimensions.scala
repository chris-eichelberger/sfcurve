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

  trait Discretizor{
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

  trait SpaceFillingCurve extends Discretizor with RangeConsolidator {
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

  }
}
