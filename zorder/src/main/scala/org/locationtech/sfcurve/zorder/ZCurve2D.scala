/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve._

import scala.util.Try

/** Represents a 2D Z order curve that we will use for benchmarking purposes in the early stages.
  *
  * Assumes, because of the dependence on the Z2 object to do the real range-query identification,
  * that query ranges will already be returned in order without any additional post-processing.
  *
  * NB:  Because of the parent contract, this will always be a square Hilbert curve.
  *
  * @param    resolution     The number of cells in each dimension of the grid space that will be indexed.
  */
class ZCurve2D(resolution: Int) extends SpaceFillingCurve2D(Dimensions.bitsFromCardinality(resolution).toInt)
  with IdentityRangeConsolidator {

  val name: String = "ZOrder"

  override def fold(subordinates: Seq[Long]): Long = Z2(subordinates.head.toInt, subordinates.last.toInt).z

  override def unfold(index: Long): Vector[Long] = {
    val (ix, iy) = Z2(index).decode
    Vector(ix, iy)
  }

  override def indexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long], hints: Option[RangeComputeHints]): Seq[IndexRange] = {
    require(lowerCorner.size == 2)
    require(upperCorner.size == 2)
    val min = Z2(lowerCorner.head.toInt, lowerCorner.last.toInt)
    val max = Z2(upperCorner.head.toInt, upperCorner.last.toInt)

    val maxRecurse = for {
      hint    <- hints
      recurse <- Option(hint.get(ZCurve2D.MAX_RECURSE))
      asInt   <- Try(recurse.asInstanceOf[Int]).toOption
    } yield {
      asInt
    }

    Z2.zranges(Array(ZRange(min, max)), maxRecurse = maxRecurse)
  }

}

object ZCurve2D {
  val DEFAULT_MAX_RECURSION = 32
  val MAX_RECURSE = "zorder.max.recurse"

  def hints(maxRecurse: Int): Option[RangeComputeHints] = {
    val hints = new RangeComputeHints()
    hints.put(MAX_RECURSE, Int.box(maxRecurse))
    Some(hints)
  }
}