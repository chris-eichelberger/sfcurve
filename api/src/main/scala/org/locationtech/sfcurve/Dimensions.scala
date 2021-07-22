package org.locationtech.sfcurve

import java.util.Date

import org.locationtech.sfcurve.Utilities.CartesianProductIterable

/**
  * This object contains definitions that are at the heart of the curve abstraction
  * and composition.  The main idea is that every curve innately consists of two
  * different types of invertible transformations:  discretization and permutation.
  *
  * Discretization is responsible for translating between the (potentially continuous)
  * user space and the (always discrete) index space.  Given a point in user space,
  * it must always be possible to identify which grid-cell it belongs in.  The
  * inverse is less clean:  Given a grid-cell, it maps back to a cell, or a collection
  * of ranges per dimension, in user space.
  *
  * Permutation is how the space-filling curve orders grid-cell indexes.
  *
  * An important implication of separating these two capabilities is that we can
  * compose (think:  nest) space-filling curves.  Instead of baking WGS84 into
  * a Z-order curve, recognize that the curve is more properly expressed as:
  *
  *   Z-order (implicitly 9 bits)
  *   |
  *   +-- latitude:  4 bits; [-90, 90)
  *   |
  *   +-- longitude:  5 bits; [-180, 180)
  *
  * In this form, the dimensions are functions nested beneath the Z-order curve,
  * itself a function.  That gives us the flexibility to nest curves like this:
  *
  *   Row-major (implicitly 20 bits)
  *   |
  *   +-- ms since epoch:  11 bits; [-MaxLong, MaxLong)
  *   |
  *   +-- Z (implicitly 9 bits)
  *       |
  *       +-- latitude:  4 bits; [-90, 90)
  *       |
  *       +-- longitude:  5 bits; [-180, 180)
  *
  * TANSTAAFL:  Computing the (inverse) index of a single point is relatively fast,
  * even in a composed curve.  The expense is incurred in computing the query ranges.
  * It's expensive because each child curve can return its own (lengthy) list of
  * index ranges, and each combination of those ranges has to be run against the parent
  * curve to identify the qualifying index ranges.  Reusing the composed-
  * curve example above -- R(t, Z(y, x)) -- think about computing the range query
  * for "the state of Virginia over May 7, 2021 to August 8, 2021"...
  *
  * 1.  At the top level, the row-major curve passes the time portion of the query
  *     into the child that handles time.  Assume that it comes back with a single
  *     range, [1998-2001].  This is fast, being a single dimension.
  * 2.  The row-major curve passes the geographic portion of the query into the
  *     Z-order curve.  Assume that it comes back with three (gapless) ranges:
  *     [156-158, 161-162, 173].  This is relatively fast for a Z-order curve.
  * 3.  The row-major curve that was originally given just two constraints -- a
  *     polygon and a time interval -- now has three separate queries to plan,
  *     the cross-product of all ranges passed back from its children.  For a row-
  *     major curve, this example isn't too bad, but if the first child had used
  *     more than one dimension, you would expect multiple ranges to come back
  *     frequently, which would mean that the parent curve has a much larger
  *     number of ranges to plan than any of its children.  This will become slow
  *     for all but the easy/simple cases (such as this example).
  *
  */
object Dimensions {
  case class Extent[T](min: T, max: T, incMin: Boolean = true, incMax: Boolean = false)(implicit ev: Ordering[T]) {
    def contains(v: Any): Boolean = v match { case value: T =>
      //val value: T = v.asInstanceOf[T]
      if (ev.lt(value, min)) return false
      if (ev.gt(value, max)) return false
      if (ev.equiv(value, min) && !incMin) return false
      if (ev.equiv(value, max) && !incMax) return false
      true
    }

    def overlaps(that: Extent[T]): Boolean = {
      // easy rejection
      if (ev.lt(that.max, min)) return false
      if (ev.gt(that.min, max)) return false

      // endpoint rejection
      if (ev.eq(that.max, min) && !(that.incMax && incMin)) return false
      if (ev.eq(that.min, max) && !(that.incMin && incMax)) return false

      // if you get this far, you overlap
      true
    }

    // ignores openness
    def isPoint: Boolean = ev.equiv(min, max)
  }

  implicit val dateOrdering: Ordering[Date] = new Ordering[Date] {
    def compare(a: Date, b: Date): Int = {
      a.getTime.compare(b.getTime)
    }
  }

  trait DimensionLike[T] {
    def toBin(value: T, extent: Extent[T], cardinality: Long): Long
    def toExtent(idx: Long, extent: Extent[T], cardinality: Long): Extent[T]
  }

  implicit object LongDimensionLike extends DimensionLike[Long] {
    def toBin(value: Long, extent: Extent[Long], cardinality: Long): Long = {
      //require(extent.contains(value))
      Math.max(0, Math.min(cardinality - 1, Math.floor((value - extent.min) * cardinality.toDouble / (extent.max - extent.min)).toLong))
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
      //require(extent.contains(value))
      Math.max(0, Math.min(cardinality - 1, ((value - extent.min) * cardinality / (extent.max - extent.min)).toLong))
      //((x + 180) * (precision - 1) / 360d).toLong
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
    def index(values: Seq[Any]): Long
    def inverseIndex(index: Long): Cell
    def queryRanges(extents: Seq[Option[Extent[_]]], hints: Option[RangeComputeHints] = None): Iterator[IndexRange]
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

    def index(values: Seq[Any]): Long = {
      require(values.size == 1)
      values.head match {
        case value: T => toBin(normalize(value))
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
    def queryRanges(extents: Seq[Option[Extent[_]]], hints: Option[RangeComputeHints] = None): Iterator[IndexRange] = {
      require(extents.size == 1)
      val subExtent = extents.head
      if (subExtent.isEmpty) return Seq(fullIndexRange).iterator
      if (subExtent.contains(extent)) return Seq(fullIndexRange).iterator
      Seq(CoveredRange(toBin(subExtent.get.min.asInstanceOf[T]), toBin(subExtent.get.max.asInstanceOf[T]))).iterator
    }
  }

  class DoubleDimension(min: Double, max: Double, val cardinality: Long) extends Dimension[Double] {
    require(max > min)
    val ev: DimensionLike[Double] = implicitly[DimensionLike[Double]]
    val extent: Extent[Double] = Extent[Double](min, max, incMin = true, incMax = false)
    val length: Double = max - min
    override def normalize(value: Double): Double = {
      var v = value
      while (v < min) v += length
      while (v >= max) v -= length
      v
    }
    override def toString: String = {
      s"DoubleDimension:  [$min, $max), cardinality $cardinality)"
    }
  }

  // defined over [-180.0, 180.0)
  case class Longitude(override val cardinality: Long) extends DoubleDimension(-180.0, 180.0, cardinality)

  // defined over [-90.0, 90.0)
  case class Latitude(override val cardinality: Long) extends DoubleDimension(-90.0, 90.0, cardinality)

  // defined over [-180.0, 180.0)
  case class NonNegativeLongitude(override val cardinality: Long) extends DoubleDimension(0.0, 360.0, cardinality)

  // defined over [-90.0, 90.0)
  case class NonNegativeLatitude(override val cardinality: Long) extends DoubleDimension(0.0, 180.0, cardinality)

  // defined over [-90.0, 90.0)
  case class AltitudeInMeters(override val cardinality: Long) extends DoubleDimension(-11000.0, 17000.0, cardinality)

  case class DateDimension(extent: Extent[Date], cardinality: Long) extends Dimension[Date] {
    val ev: DimensionLike[Date] = implicitly[DimensionLike[Date]]
    override def normalize(value: Date): Date = {
      if (!extent.contains(value)) throw new Exception(s"Invalid date within the date dimension:  $value outside $extent")
      value
    }
  }
  def DateDimension(minDate: Date, maxDate: Date, cardinality: Long): DateDimension =
    DateDimension(Extent[Date](minDate, maxDate), cardinality)

  // defined over [1970-01-01, 2070-01-01)
  def ExampleEra(cardinality: Long): DateDimension = {
    val d0 = new Date(70, 0, 1)
    val d1 = new Date(170, 0, 1)
    DateDimension(d0, d1, cardinality)
  }

  case class Cell(extents: Vector[Extent[_]] = Vector()) {
    val numDimensions: Int = extents.size

    def contains(point: Vector[_]): Boolean = {
      if (point.size != numDimensions) return false
      extents.zip(point).forall {
        case (extent, coord) => extent.contains(coord)
      }
    }

    def +(that: Cell): Cell = Cell(extents ++ that.extents)
  }

  implicit def hintsToGap(hintsOpt: Option[RangeComputeHints]): Long = {
    if (hintsOpt == null) return 0
    val hint = hintsOpt.map(hints => hints.get(GapMergedIndexRange.HintsKeyMapGap)).getOrElse(0)
    if (hint == null) return 0
    hint.toString.toLong
  }

  def bitsFromCardinality(cardinality: Long): Long =
    Math.ceil(Math.log(cardinality.toDouble) / Math.log(2.0)).toLong

  trait SpaceFillingCurve extends Discretizor with RangeConsolidator {
    def children: Vector[Discretizor]

    //require(children != null && children.nonEmpty)

    // how many possible values can the `index` method return?
    def cardinality: Long = children.map(_.cardinality).product

    // how many dimension leaf-nodes are there?
    def arity: Int = children.map(_.arity).sum

    // fetch the leaf-node discretizors (dimensions, typically)
    def leaves: Vector[Discretizor] = children.map {
      case _: SpaceFillingCurve => leaves
      case dimension: Dimension[_] => Vector(dimension)
      case other => throw new Exception(s"Invalid SFC child type:  ${other.getClass}")
    }.reduce(_ ++ _)

    // these routines are the heart of the SFC;
    // they will be made concrete per descendant curve
    def fold(subordinates: Seq[Long]): Long
    def unfold(index: Long): Vector[Long]

    // ideally, this will emit qualifying index ranges in order, but that may not always
    // be possible...
    def indexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long], hints: Option[RangeComputeHints] = None): Seq[IndexRange]

    def index(values: Seq[Any]): Long = {
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

    // queries ONLY a SINGLE rectangle, converting user-space ranges to index-space ranges recursively
    // NB:  "maxGap" is passed through the (optional) RangeComputeHints
    def queryRanges(extents: Seq[Option[Extent[_]]], hints: Option[RangeComputeHints] = None): Iterator[IndexRange] = {
      // assume that unqueried dimensions are either:
      // 1.  an option set to its full extent; or
      // 2.  a None
      require(extents.size == arity)

      // extract the maximum allowable gap between ranges (defaults to zero)
      val maxGap: Long = hints

      // convert the per-dimension extents into (flat) index ranges per child
      // (defer to your children for this work, be they SFCs or Dimensions)
      val childIndexes: Seq[Seq[IndexRange]] = children.foldLeft((extents, Seq[Seq[IndexRange]]()))((acc, child) => acc match {
        case (extentsLeft: Seq[Option[Extent[_]]], subsSoFar: Seq[IndexRange]) =>
          val childExtents: Seq[Option[Extent[_]]] = extentsLeft.take(child.arity)
          val childRanges: Seq[IndexRange] = child.queryRanges(childExtents, hints).toSeq
          (extentsLeft.drop(child.arity), subsSoFar :+ childRanges)
      })._2

      // it's possible that a single child's conversion of its user ranges returns multiple query ranges --
      // think of R(Z(x, y), t) when presented with ([-10, 10], [0, 5], [now, now-5min]); the
      // Z curve is discontinuous across the x- and y- box, returning at least two separate ranges --
      // so plan across all combinations of child index ranges
      val rectangleIterator: Iterator[Seq[IndexRange]] =
        CartesianProductIterable(childIndexes).iterator.map(_.asInstanceOf[Seq[IndexRange]])
      val rangesToConsolidate: Iterator[IndexRange] = rectangleIterator.flatMap(indexRangeSeq => {
        val lowerCorner = indexRangeSeq.map(_.lower)
        val upperCorner = indexRangeSeq.map(_.upper)
        indexRanges(lowerCorner, upperCorner, hints)
      })

      // roll up the child index-ranges into a single index-range
      // TODO:  short-cut if there was only one rectangle; these ranges will already be in order
      consolidateRanges(rangesToConsolidate, maxGap)
    }

  }
}
