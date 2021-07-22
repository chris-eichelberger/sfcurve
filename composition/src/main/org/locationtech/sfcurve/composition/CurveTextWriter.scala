package org.locationtech.sfcurve.composition

import org.locationtech.sfcurve.Dimensions.{Dimension, Discretizor, SpaceFillingCurve}
import org.locationtech.sfcurve.Utilities.CartesianProductIterable

/**
 * This is a utility class that writes a (composed) space-filling curve in
 * plain-text form to stdout.  This is useful for debugging and for some
 * simple visualization.
 */
object CurveTextWriter {
  type Indexer = Seq[Long] => Long
  type OutputSink = String => Unit

  // draws a single, 1D or 2D slice of the index-space
  private def writeGrid(write: OutputSink, higherCoord: Seq[Long], width: Int, height: Int, indexer: Indexer, cellWidth: Int = 5): Unit = {
    // if there are higher coordinates you're iterating over, show those first for this slice
    if (higherCoord.nonEmpty) {
      val numHigherDims: Int = higherCoord.size
      higherCoord.foldLeft(numHigherDims)((idx, coord) => {
        write(s"Dimension #$idx:  $coord")
        idx - 1
      })
    }

    write(s"Slice size:  $width columns by $height rows")

    // draw this grid slice
    val HR: String = "+" + ("-" * (2 + cellWidth) + "+") * width.toInt
    for (row <- (0 until height).reverse) {
      write(HR)
      write((0 until width).map(col => ("%" + cellWidth.toString + "d").format(indexer(Seq[Long](col, row) ++ higherCoord))).mkString("| ", " | ", " |"))
    }
    write(HR)
  }

  /**
   * This HOF is needed because the SFC contract is only meant to recurse for raw values,
   * not for discretized indexes.  That is, the SFC's `fold` expects a list that contains
   * one Long per *child*, whereas `index` expects a list that contains one value per *dimension*.
   *
   * This is an artificial case that should only be relevant inside unit tests, which is why
   * this function isn't simply included into the SFC class itself.
   *
   * @param discretizor the curve or dimension
   * @param coord the list of per-dimension discretized values
   * @return the `fold`ed result of using these per-dimension values up through the curve
   */
  def index(discretizor: Discretizor)(coord: Seq[Long]): Long = {
    discretizor match {
      case curve: SpaceFillingCurve =>
        val toFold: Seq[Long] = curve.children.foldLeft(coord, Seq[Long]())((acc, child) => acc match {
          case (coordRemainder, toFoldSoFar) =>
            (
              coordRemainder.drop(child.arity),
              toFoldSoFar :+ index(child)(coordRemainder)
            )
        })._2
        curve.fold(toFold)
      case dim: Dimension[_] =>
        coord.head
    }
  }

  // the first dimension iterates fastest; the last dimension iterates slowest
  def writeText(curve: SpaceFillingCurve, write: OutputSink = println): Unit = {
    val leaves: Vector[Discretizor] = curve.leaves

    write(s"Writing text representation of curve:  $curve")
    write(s"Leaves:  ${leaves.map(_.toString).mkString("; ")}")
    write(s"Cardinality:  ${curve.cardinality}")

    // some curves are just too big to print
    require(curve.cardinality <= 8192)
    require(leaves.forall(_.cardinality <= 64))

    require(leaves.nonEmpty)
    val width: Int = leaves(0).cardinality.toInt
    val height: Int = if (leaves.length > 1) leaves(1).cardinality.toInt else 1

    val cellWidth: Int = Math.ceil(Math.log10(curve.cardinality)).toInt

    val higherCardinalities: Vector[Vector[Long]] = leaves.drop(2).map(child => (0 until child.cardinality.toInt).map(_.toLong).toVector)
    write(s"Higher cardinalities:  ${higherCardinalities.map(_.toString).mkString("[", ", ", "]")}")
    CartesianProductIterable(higherCardinalities).iterator.foreach { higherCoord =>
      writeGrid(write, higherCoord.asInstanceOf[Seq[Long]], width, height, index(curve), cellWidth)
    }
  }
}
