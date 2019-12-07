package org.locationtech.sfcurve.zorder

import java.io.{BufferedWriter, FileWriter, OutputStreamWriter, PrintWriter}

import org.locationtech.sfcurve.Utilities.CartesianProductIterable

import scala.collection.mutable.{Map => MutableMap, Queue => MutableQueue, Set => MutableSet, Stack => MutableStack}

object Gray extends App {

  case class Ring(cardinality: Int) {
    def adjacents(i: Int): Set[Int] = Set(
      Math.abs((i - 1) % cardinality),
      Math.abs((i + 1) % cardinality)
    ) - i

    def startingSeq: Seq[Int] = {
      val start = 0
      val adjs = adjacents(start)
      Set(start).union(adjs).toList.sorted
    }

    def valid(a: Int): Boolean = (a >= 0 && a < cardinality)

    def areNotFar(a: Int, b: Int): Boolean = valid(a) && valid(b) &&
      (Math.abs(a - b) < 2 || (Math.min(a, b) == 0 && Math.max(a, b) == cardinality - 1))

    def areDifferentAndAdjacent(a: Int, b: Int): Boolean = (a != b) && areNotFar(a, b)

    def seq: Seq[Int] = (0 until cardinality).toSeq
  }

  case class Code(values: Vector[Int])

  case class CodeScheme(rings: Vector[Ring]) {
    def codes: Seq[Code] = {
      val places: Seq[Seq[Int]] = rings.map(_.seq)
      CartesianProductIterable(places).toSeq.map {
        case s: Seq[Int] => Code(s.toVector)
      }
    }

    def adjacents(code: Code): Set[Code] = {
      val roots: Seq[Seq[Int]] = rings.zip(code.values).map {
        case (ring, value) => ring.adjacents(value).toSeq
      }
      CartesianProductIterable(roots).map {
        case s: Seq[Int] => Code(s.toVector)
        case _ => throw new Exception("Unexpected type")
      }.toSet
    }

    def areDifferentAndAdjacent(a: Code, b: Code): Boolean = {
      rings.zip(a.values).zip(b.values).map {
        case ((ring, aValue), bValue) =>
          if (ring.areNotFar(aValue, bValue)) {
            if (aValue == bValue) 0
            else 1
          } else 2
      }.sum == 1
    }
  }

  case class Node(id: Int) {
    def toDot(encodingOpt: Option[Map[Node, Code]]): String = {
      val label = if (encodingOpt.isDefined) {
        val codeOpt: Option[Code] = encodingOpt.get.get(this)
        if (codeOpt.isEmpty) throw new Exception(s"Undefined code for node $id")
        codeOpt.get.values.mkString(",")
      } else ""
      s"""node_$id [ label=\"$label\" shape=\"circle\" width=0.3 ]"""
    }
  }

  case class Edge(fromNode: Node, toNode: Node) {
    require(fromNode.id < toNode.id)
    def toDot: String = s"node_${fromNode.id} -- node_${toNode.id}"
  }

  case class Graph(nodes: Vector[Node], edges: Vector[Edge]) {
    def size: Int = nodes.size

    def connections(node: Node): Set[Node] =
      edges.filter(e => e.fromNode == node || e.toNode == node).map(e => if (e.fromNode == node) e.toNode else e.fromNode).toSet

    val connectionDistances: MutableMap[Node, Map[Node, Int]] = MutableMap.empty[Node, Map[Node, Int]]

    def getConnectionDistances(start: Node): Map[Node, Int] = {
      if (connectionDistances.contains(start)) return connectionDistances(start)

      var cnxs: Set[Node] = connections(start)
      if (cnxs.size < 1) Map.empty[Node, Int]

      val minimums = MutableMap.empty[Node, Int]
      minimums.put(start, 0)

      def acc(a: Node, distance: Int): Unit = {
        if (minimums.contains(a)) {
          minimums.put(a, Math.min(minimums(a), distance))
        } else {
          minimums.put(a, distance)
        }
      }

      val visited: MutableSet[Node] = MutableSet[Node]()
      val queue: MutableQueue[(Node, Int)] = MutableQueue.empty[(Node, Int)]
      connections(start).foreach(other => queue.enqueue((other, 1)))
      while (queue.nonEmpty) {
        val (node, distance) = queue.dequeue()
        acc(node, distance)
        visited.add(node)
        cnxs = connections(node)
        cnxs.diff(visited).diff(queue.map(_._1).toSet).foreach(n => queue.enqueue((n, distance + 1)))
      }

      //println(s"Cx($start):  $minimums")

      connectionDistances.put(start, minimums.toMap)

      minimums.toMap
    }

    def getConnectionsAtDistance(start: Node, d: Int): Set[Node] = {
      getConnectionDistances(start).filter { case (k, v) => v == d }.keySet
    }

    def connectedness(node: Node): Double = {
      val base = 1.0 / (nodes.length.toDouble - 1.0)
      (1 until nodes.length).foldLeft(0.0)((acc, distance) => {
        val count = getConnectionsAtDistance(node, distance).size.toDouble
        val weight = Math.pow(base, distance.toDouble)
        acc + (count * weight)
      })
    }

    def toDot(encodingOpt: Option[Map[Node, Code]] = None): String =
      "graph G {\n" +
      nodes.map(node => node.toDot(encodingOpt) + f" cnx=${connectedness(node)}%1.4f").mkString("  ", "\n  ", "\n") +
      edges.map(_.toDot).mkString("  ", "\n  ", "\n") +
      "}"

    def firstNode: Node = nodes.head

    // assumes that the graph is connected
    def getEncoding: Option[Map[Node, Code]] = {
      val cardinalities: Seq[Seq[Int]] = for (i <- 1 to nodes.size) yield Seq(3, 2)
      val ringSizes: Seq[Seq[_]] = CartesianProductIterable(cardinalities).toSeq
      val schemas: Seq[CodeScheme] = ringSizes.map {
        case s: Seq[Int] =>
          val rings = s.map(cardinality => Ring(cardinality)).toVector
          CodeScheme(rings)
      }

      def findEncoding(scheme: CodeScheme): Option[Map[Node, Code]] = {
        val NoMapping = Map.empty[Node,Code]
        //println(s"    findEncoding $scheme")
        val codes = scheme.codes.toList
        //println(s"      codes:  $codes")

        def solve(nodesLeft: Set[Node], codesLeft: Set[Code], mapSoFar: Map[Node, Code]): Option[Map[Node, Code]] = {
          if (nodesLeft.isEmpty) return Option(mapSoFar)

          def acceptable(node: Node, code: Code): Boolean = {
            connections(node).forall(otherNode => mapSoFar.get(otherNode) match {
              case Some(otherCode) =>
                scheme.areDifferentAndAdjacent(otherCode, code)
              case None =>
                true
            })
          }

          val seqs: Seq[Seq[_]] = Seq(nodesLeft.toSeq, codesLeft.toSeq)
          val choices: Iterator[Seq[_]] = CartesianProductIterable(seqs).iterator

          choices.flatMap {
            case Seq(node: Node, code: Code) =>
              if (acceptable(node, code)) {
                //println(s"        ACCEPTABLE:  $node -> $code")
                val candidate: Option[Map[Node, Code]] = solve(nodesLeft - node, codesLeft - code, mapSoFar + (node -> code))
                candidate.filter(_.size == nodes.size)
              }
              else {
                //println(s"        unacceptable:  $node -> $code")
                None
              }
          }.find(_.size == nodes.size)
        }

        val startNode: Node = nodes.head
        val startCode: Code = codes.head
        val startMap: Map[Node, Code] = Map[Node, Code](startNode -> startCode)

        solve(nodes.toSet - startNode, codes.toSet - startCode, startMap)
      }

      schemas.map(scheme => findEncoding(scheme)).find(solnOpt => solnOpt.isDefined && solnOpt.get.size == nodes.size).flatten
    }

    def isConnected: Boolean = {
      var node = firstNode
      var cnxs: Set[Node] = connections(node)
      if (cnxs.size < 1) return false
      val visited: MutableSet[Node] = MutableSet[Node]()
      val stack: MutableStack[Node] = MutableStack[Node](node)

      while (stack.nonEmpty) {
        node = stack.pop()
        visited.add(node)
        cnxs = connections(node)
        cnxs.diff(visited).diff(stack.toSet).foreach(stack.push)
      }

      visited.size == nodes.size
    }

    def canonicalForm: Graph = {
      def nodeComparator(a: Node, b: Node): Boolean = {
        val cA = connectedness(a)
        val cB = connectedness(b)
        if (cA > cB) return false
        if (cA < cB) return true
        require(cA == cB)
        true
      }
      def edgeComparator(a: Edge, b: Edge): Boolean = {
        if (a.fromNode.id > b.fromNode.id) return false
        if (a.fromNode.id < b.fromNode.id) return true
        require(a.fromNode.id == b.fromNode.id)
        if (a.toNode.id > b.toNode.id) return false
        if (a.toNode.id < b.toNode.id) return true
        require(a.toNode.id == b.toNode.id)
        throw new Exception(s"Duplicate edges:  $a == $b")
      }
      val sortedNodeIds = nodes.toList.sortWith(nodeComparator).map(_.id)
      //sortedNodeIds.foreach(id => println(f"  node ID $id%s, cnx ${connectedness(nodes.find(_.id == id).get)}%1.4f"))
      val origNodeIds = nodes.indices.toList
      val remapIds: Map[Int, Int] = sortedNodeIds.zip(origNodeIds).toMap
      val newNodes: Vector[Node] = origNodeIds.map(i => Node(i)).toVector
      val newEdges: Vector[Edge] = edges.map(edge => {
        val a = remapIds(edge.fromNode.id)
        val b = remapIds(edge.toNode.id)
        Edge(
          newNodes(Math.min(a, b)),
          newNodes(Math.max(a, b))
        )
      }).sortWith(edgeComparator)
      Graph(newNodes, newEdges)
    }

    def render(outf: String, encodingOpt: Option[Map[Node, Code]]): Unit = {
      var pw: PrintWriter = null
      try {
        pw = new PrintWriter(new FileWriter("testCombination.dot"))
        pw.println(toDot(encodingOpt))
        pw.close()

        // render the output graph
        val cmd: String = s"c:\\Program Files (x86)\\Graphviz2.38\\bin\\neato.exe -T png testCombination.dot -o $outf"
        val p: Process = Runtime.getRuntime.exec(cmd)
        require(p.waitFor() == 0)
      } finally {
        pw.close()
      }
    }
  }

  object Graph {
    def apply(n: Long, numBits: Int, combinations: Long, bits: Seq[Int]): Graph = {
      val nodes: Vector[Node] = (0 until n.toInt).map(i => Node(i)).toVector
      var fromIndex = 0
      var toIndex = 1
      val edges: Vector[Edge] = (0 until numBits).foldLeft(Vector[Edge]())((acc, j) => {
        val next = if (bits(j) > 0) acc :+ Edge(nodes(fromIndex), nodes(toIndex)) else acc
        toIndex += 1
        if (toIndex >= n) {
          fromIndex += 1
          toIndex = fromIndex + 1
        }
        next
      })
      Graph(nodes, edges)
    }

    def simplifyEncoding(encodingOpt: Option[Map[Node, Code]]): Option[Map[Node, Code]] = {
      encodingOpt.map(encoding => {
        val codes: Seq[Code] = encoding.values.toSeq
        val numPlaces = codes.head.values.size
        val placesToKeep: Vector[Boolean] = (0 until numPlaces).map(place => {
          val values: Int = codes.map(_.values(place)).toSet.size
          values > 1
        }).toVector
        encoding.map {
          case (node, code) =>
            val newValues: Vector[Int] = code.values.zipWithIndex.flatMap {
              case (value, place) => if (placesToKeep(place)) Option(value) else None
            }
            node -> Code(newValues)
        }
      })
    }
  }

  def generate(n: Long): Long = {
    var numConnected: Long = 0

    val numBits: Int = (((n - 1) * n) >> 1).toInt
    println(s"Total possible configuration bits:  $numBits")

    val combinations: Long = 1L << numBits
    println(s"Total number of possible configurations:  $combinations")

    val ucombs = MutableSet[String]()

    for (i <- 0L until combinations) {
      val bits: Seq[Int] = i.toBinaryString.split("").reverse.padTo(numBits, "0").reverseMap(_.toInt)
      val graph = Graph(n, numBits, combinations, bits)
      val isConn = graph.isConnected
      println(s"  ${n}-node combination $i / $combinations:  $bits, connected $isConn")

      if (isConn) {
        val encodingOpt = Graph.simplifyEncoding(graph.getEncoding)
        val comb = graph.canonicalForm.toDot().replaceAll("\n", "  ")
        if (!ucombs.contains(comb)) {
          if (encodingOpt.isDefined) println(s"    encoding:  ${encodingOpt.get}")
          else println(s"    NO ENCODING FOUND")
          graph.render(s"testCombination_n${n}_c${i}.png", encodingOpt)
        }
        ucombs.add(comb)
      }
    }

    numConnected
  }

  (1 to 5).foreach(n => generate(n))
}
