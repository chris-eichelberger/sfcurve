package org.locationtech.sfcurve

import org.scalatest.{FunSpec, Matchers}

import GapMergedIndexRange._

class RangeConsolidatorSpec extends FunSpec with Matchers {
  implicit def pairToRange(pair: (Long, Long)): IndexRange = IndexRange(pair._1, pair._2, contained = false)

  val rawRangesOne: Seq[(Long, Long)] = Seq(
    (3, 5),
    (8, 8),
    (0, 1),
    (4, 6)
  )

  val rawRangesTwo: Seq[(Long, Long)] = Seq(
    (10, 100),
    (200, 205),
    (300, 390),
    (400, 450),
    (500, 500),
    (700, 750),
    (800, 815),
    (900, 1000)
  )

  val rawRangesThree: Seq[(Long, Long)] = Seq(
    (500, 600),
    (545, 555),
    (450, 550),
    (550, 650)
  )

  val rawRangesFour: Seq[(Long, Long)] = (16L to 31L).map(i => (i, i))

  def consolidateRanges(ranges: Seq[(Long, Long)], name: String, maxGap: Long): Int = {
    val netRanges = ranges.foldLeft(IndexRangeTree())((acc, range) =>
      acc.add(GapMergedIndexRange(range.lower, range.upper, maxGap))).
      inOrderExtentIterator.toSeq

    println(s"\nNet ranges '$name' with gap $maxGap...")
    netRanges.foreach(println)

    netRanges.size
  }

  def consolidateRangesBothWays(ranges: Seq[(Long, Long)], name: String, maxGap: Long): Int = {
    val forwards = consolidateRanges(ranges, name, maxGap)
    val backwards = consolidateRanges(ranges.reverse, s"$name, Reverse", maxGap)
    if (forwards == backwards) forwards else -1
  }

  describe("in-memory range consolidator") {
    it("should merge the ranges with 0 gap") {
      consolidateRangesBothWays(rawRangesOne, "One", 0) should equal(3)
      consolidateRangesBothWays(rawRangesTwo, "Two", 0) should equal(8)
      consolidateRangesBothWays(rawRangesThree, "Three", 0) should equal(1)
      consolidateRangesBothWays(rawRangesFour, "Four", 0) should equal(1)
    }

    it("should merge the ranges with 50 gap") {
      consolidateRangesBothWays(rawRangesOne, "One", 50) should equal(1)
      consolidateRangesBothWays(rawRangesTwo, "Two", 50) should equal(5)
      consolidateRangesBothWays(rawRangesThree, "Three", 50) should equal(1)
      consolidateRangesBothWays(rawRangesFour, "Four", 50) should equal(1)
    }

    it("should merge the ranges with 150 gap") {
      consolidateRangesBothWays(rawRangesOne, "One", 150) should equal(1)
      consolidateRangesBothWays(rawRangesTwo, "Two", 150) should equal(2)
      consolidateRangesBothWays(rawRangesThree, "Three", 150) should equal(1)
      consolidateRangesBothWays(rawRangesFour, "Four", 150) should equal(1)
    }
  }
}
