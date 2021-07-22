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

/** Represents a 3D Z order curve that we will use for benchmarking purposes in the early stages.
  *
  * Assumes, because of the dependence on the Z3 object to do the real range-query identification,
  * that query ranges will already be returned in order without any additional post-processing.
  *
  * @param    resolution     The number of cells in each dimension of the grid space that will be indexed.
  */
class ZCurve3D(resolution: Int) extends SpaceFillingCurve3D(Dimensions.bitsFromCardinality(resolution).toInt)
  with IdentityRangeConsolidator {

  override def fold(subordinates: Seq[Long]): Long = Z3(subordinates.head.toInt, subordinates(1).toInt, subordinates.last.toInt).z

  override def unfold(index: Long): Vector[Long] = {
    val (ix, iy, iz) = Z3(index).decode
    Vector(ix, iy, iz)
  }

  override def indexRanges(lowerCorner: Seq[Long], upperCorner: Seq[Long], hints: Option[RangeComputeHints]): Seq[IndexRange] = {
    require(lowerCorner.size == 3)
    require(upperCorner.size == 3)
    val min = Z3(lowerCorner.head.toInt, lowerCorner(1).toInt, lowerCorner.last.toInt)
    val max = Z3(upperCorner.head.toInt, upperCorner(1).toInt, upperCorner.last.toInt)

    val maxRecurse = for {
      hint    <- hints
      recurse <- Option(hint.get(ZCurve3D.MAX_RECURSE))
      asInt   <- Try(recurse.asInstanceOf[Int]).toOption
    } yield {
      asInt
    }

    Z3.zranges(Array(ZRange(min, max)), maxRecurse = maxRecurse)
  }

}

object ZCurve3D {
  val DEFAULT_MAX_RECURSION = 32
  val MAX_RECURSE = "zorder.max.recurse"

  def hints(maxRecurse: Int): Option[RangeComputeHints] = {
    val hints = new RangeComputeHints()
    hints.put(MAX_RECURSE, Int.box(maxRecurse))
    Some(hints)
  }
}