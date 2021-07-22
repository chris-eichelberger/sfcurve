package org.locationtech.sfcurve

import org.locationtech.sfcurve.Dimensions.{Cell, Dimension, Discretizor, SpaceFillingCurve}
import org.locationtech.sfcurve.Utilities.CartesianProductIterable

import scala.Console.err

/**
 * This object is used to describe tests that should ONLY BE RUN ON SMALL CURVES
 * because it does exhaustive testing about the indexing to make sure tha the
 * curve represents a valid permutation.
 */
case class CurveValidation(curve: SpaceFillingCurve) {
  // the curve is said to be "pure" if all children are dimension discretizors;
  // contrast with "composed"
  def isPure: Boolean = curve.children.forall {
    case subcurve: SpaceFillingCurve => false
    case dimension: Dimension[_] => true
  }

  def isSizeOkay: Boolean = {
   curve.cardinality <= (1 << 20) &&
     curve.children.nonEmpty
  }

  def isValidPermutation[T: Ordering](set: Set[T], min: T, max: T): Boolean = {
    if (set == null) {
      err.println("set was `null`")
      return false
    }

    if (set.size != curve.cardinality) {
      err.println(s"Set contained ${set.size} element(s) when ${curve.cardinality} were expected")
      return false
    }

    if (set.min != min) {
      err.println(s"Set contained minimum ${set.min} when $min was expected")
      return false
    }

    if (set.max != max) {
      err.println(s"Set contained maximum ${set.max} when $max was expected")
      return false
    }

    true
  }

  // given all valid combinations of child indexes, does the top-level
  // curve produce a valid permutation?
  def topFoldWorks: Boolean = {
    val cardinalities: Vector[Vector[Long]] = curve.children.map(child => (0 until child.cardinality.toInt).map(_.toLong).toVector)
    val set: Set[Long] =
      CartesianProductIterable(cardinalities).iterator.map { coord =>
        curve.fold(coord.asInstanceOf[Seq[Long]])
      }.toSet

    isValidPermutation[Long](set, 0, curve.cardinality - 1L)
  }

  // given all valid combinations of output indexes, does the top-level
  // curve produce a valid permutation from the inverses?
  def topUnfoldWorks: Boolean = {
    val set: Set[Seq[Long]] = (0 until curve.cardinality.toInt).map(i => curve.unfold(i.toLong)).toSet

    isValidPermutation[Seq[Long]](set, Seq(0, 0, 0), curve.children.map(_.cardinality - 1L))(CurveValidation.CoordinateOrdering)
  }

  // this tests the full curve, pure or composed, to make sure that
  // the (inverse)indexing works as intended, at least generating a
  // valid permutation
  def roundTripsWork: Boolean = {
    var index = 0
    val mset = scala.collection.mutable.Set[Cell]()
    while (index < curve.cardinality) {
      // confirm the cell is unique
      val cell = curve.inverseIndex(index)
      if (mset.contains(cell)) {
        err.println(s"Duplicate cell encountered unexpectedly:  index $index -> cell $cell")
        return false
      }
      mset.add(cell)

      // confirm that you get the index back the other direction
      val point: Vector[Any] = cell.extents.map(_.min)
      require(cell.contains(point))
      val reIndex: Long = curve.index(point)
      if (reIndex != index) {
        err.println(s"Failed round-trip index:  index $index -> cell $cell -> point $point -> index $reIndex")
        return false
      }

      index += 1
    }

    // it *shouldn't* be possible for this to be false, right?
    mset.size == curve.cardinality
  }

  // this is the method we intend to call most often:
  // do all of the component tests pass to define a valid curve
  def isValid: Boolean = {
    isSizeOkay &&
      topFoldWorks &&
      topUnfoldWorks &&
      roundTripsWork
  }
}

object CurveValidation {
  val LT: Int = -1
  val EQ: Int = 0
  val GT: Int = 1

  object CoordinateOrdering extends Ordering[Seq[Long]] {
    override def compare(x: Seq[Long], y: Seq[Long]): Int = {
      if ((x == null) && (y == null)) return EQ
      if ((x == null) && (y != null)) return LT
      if ((x != null) && (y == null)) return GT

      if (x.size != y.size) return x.size.compareTo(y.size)

      x.zip(y).foreach {
        case (xx, yy) =>
          val comparison = xx.compareTo(yy)
          if (comparison != EQ) return comparison
      }

      EQ
    }
  }
}