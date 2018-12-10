package org.locationtech.sfcurve

object Utilities {
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
}
