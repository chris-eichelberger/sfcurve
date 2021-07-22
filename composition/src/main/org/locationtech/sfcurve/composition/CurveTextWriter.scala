package org.locationtech.sfcurve.composition

import org.locationtech.sfcurve.Dimensions.{Discretizor, SpaceFillingCurve}
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
    val HR: String = "+" + ("-" * cellWidth + "+") * width.toInt
    for (row <- 0 until height) {
      write(HR)
      "|" + ((0 until width).map(col => ("%" + cellWidth.toString + "d").format(indexer(Seq[Long](row, col) ++ higherCoord)) + "|"))
    }
    write(HR)
  }

  // the first dimension iterates fastest; the last dimension iterates slowest
  def writeText(curve: SpaceFillingCurve, write: OutputSink = println): Unit = {
    val leaves: Vector[Discretizor] = curve.leaves

    write(s"Writing text representation of curve:  $curve")
    write(s"Leaves:  ${leaves.map(_.toString).mkString("; ")}")

    // some curves are just too big to print
    require(curve.cardinality <= 8192)
    require(leaves.forall(_.cardinality <= 64))

    require(leaves.nonEmpty)
    val width: Int = leaves(0).cardinality.toInt
    val height: Int = if (leaves.length > 1) leaves(1).cardinality.toInt else 1

    val higherCardinalities: Vector[Vector[Long]] = curve.children.drop(2).map(child => (0 until child.cardinality.toInt).map(_.toLong).toVector)
    CartesianProductIterable(higherCardinalities).iterator.foreach { higherCoord =>
      writeGrid(write, higherCoord.asInstanceOf[Seq[Long]], width, height, curve.fold)
    }
  }
}
