/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2026                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#ifndef VCG_STREAM_ORDER_H
#define VCG_STREAM_ORDER_H

#include <vcg/complex/complex.h>
#include <vcg/simplex/edge/topology.h>

namespace vcg {
namespace tri {

/**
 * @brief Stream orders of a tree-shaped edge mesh.
 *
 * A stream order ranks the branches of a tree by their position in the hierarchy. The
 * two classic ones come from hydrology, but the same question arises for any branching
 * structure: vessel trees, plant and coral curve skeletons, neuron arbors.
 *
 * - **Strahler** numbers grow from the tips towards the root: a tip is 1, and a vertex
 *   takes the largest number among its children, plus one when two or more children
 *   share that largest number. A branch's number measures how deep the subtree it
 *   collects is.
 * - **Hack** numbers grow from the root towards the tips: the main stem is 1, a branch
 *   leaving it is 2, a branch leaving that one is 3, and so on. Everything depends on
 *   which child continues the stem at each junction; see MainStemRule.
 *
 * The mesh is an edge mesh (vertices and edges, faces ignored) and must be a tree on the
 * component that contains the root: Root() throws on a loop. Vertices in other components
 * are not reached and get 0. Degree-2 vertices simply pass their number along, so the
 * orders can be computed directly on a densely sampled curve skeleton.
 *
 * Results are per-vertex, indexed by tri::Index(m, v), sized m.vert.size(), with 0 for
 * vertices that are deleted or not reached. Requires VE adjacency (vertex and edge).
 */
template <class MeshType>
class StreamOrder
{
public:
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::EdgeType       EdgeType;
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::CoordType      CoordType;

  /// Which child continues the main stem at a junction, for Hack().
  enum MainStemRule
  {
    /// The child whose first edge deviates least from the direction of the stem so far,
    /// taken as the chord from the vertex where the current stem started to the junction.
    /// Suited to organisms that keep a growth direction, where the leading axis is the
    /// straightest one (e.g. Acropora corals).
    MinDeviationAngle,
    /// The child leading to the longest path (sum of edge lengths) down to a tip.
    /// The classical definition of Hack's order, where the main stream is the longest one.
    LongestPath
  };

  /// The tree seen from a root: parent of every vertex (-1 for the root and for vertices
  /// not reached) and the reached vertices in breadth-first order, root first.
  struct RootedTree
  {
    std::vector<int> parent;
    std::vector<int> bfs;
  };

  /// Build the rooted view of the component containing @p root. Throws
  /// MissingPreconditionException if that component is not a tree.
  static RootedTree Root(MeshType &m, VertexPointer root)
  {
    RequireVEAdjacency(m);
    if (root == nullptr || root->IsD())
      throw vcg::MissingPreconditionException("StreamOrder: the root is not a valid vertex.");

    const int rootIndex = int(tri::Index(m, root));
    RootedTree t;
    t.parent.assign(m.vert.size(), -1);
    std::vector<bool> reached(m.vert.size(), false);
    reached[rootIndex] = true;
    t.bfs.push_back(rootIndex);

    // Every edge of the component is seen twice, once from each endpoint; a tree has
    // exactly one edge fewer than it has vertices, whatever loops or parallel edges
    // would otherwise hide.
    size_t edgeEnds = 0;
    std::vector<VertexPointer> star;
    for (size_t i = 0; i < t.bfs.size(); ++i)
    {
      VertexPointer v = &m.vert[t.bfs[i]];
      edge::VVStarVE(v, star);
      edgeEnds += star.size();
      for (VertexPointer w : star)
      {
        const int wi = int(tri::Index(m, w));
        if (reached[wi]) continue;
        reached[wi] = true;
        t.parent[wi] = t.bfs[i];
        t.bfs.push_back(wi);
      }
    }
    if (edgeEnds != 2 * (t.bfs.size() - 1))
      throw vcg::MissingPreconditionException("StreamOrder: the edge mesh contains a loop, so it is not a tree.");
    return t;
  }

  /// Strahler number of every vertex.
  static std::vector<int> Strahler(MeshType &m, VertexPointer root)
  {
    const RootedTree t = Root(m, root);
    std::vector<int> order(m.vert.size(), 0);
    std::vector<int> childCount(m.vert.size(), 0);
    std::vector<int> childMax(m.vert.size(), 0);    // largest number among the children
    std::vector<int> childMaxCount(m.vert.size(), 0); // how many children reach it

    // Children before parents: every child is final when its parent is computed, so
    // the result does not depend on the order in which children are stored.
    for (auto it = t.bfs.rbegin(); it != t.bfs.rend(); ++it)
    {
      const int v = *it;
      if (childCount[v] == 0)
        order[v] = 1;
      else
        order[v] = childMax[v] + (childMaxCount[v] >= 2 ? 1 : 0);

      const int p = t.parent[v];
      if (p < 0) continue;
      ++childCount[p];
      if (order[v] > childMax[p]) { childMax[p] = order[v]; childMaxCount[p] = 1; }
      else if (order[v] == childMax[p]) ++childMaxCount[p];
    }
    return order;
  }

  /// Hack number of every vertex, with @p rule choosing the main stem at each junction.
  static std::vector<int> Hack(MeshType &m, VertexPointer root, MainStemRule rule)
  {
    const RootedTree t = Root(m, root);
    const size_t n = m.vert.size();

    std::vector<std::vector<int>> children(n);
    for (size_t i = 1; i < t.bfs.size(); ++i)
      children[t.parent[t.bfs[i]]].push_back(t.bfs[i]);

    // Length of the longest path from each vertex down to a tip. LongestPath uses it
    // everywhere; MinDeviationAngle falls back on it where there is no direction yet.
    std::vector<ScalarType> downLength(n, 0);
    for (auto it = t.bfs.rbegin(); it != t.bfs.rend(); ++it)
    {
      const int v = *it;
      for (int c : children[v])
        downLength[v] = std::max(downLength[v], downLength[c] + Distance(m.vert[v].cP(), m.vert[c].cP()));
    }
    const auto longestChild = [&](int v) {
      int best = children[v].front();
      ScalarType bestLength = -1;
      for (int c : children[v])
      {
        const ScalarType l = downLength[c] + Distance(m.vert[v].cP(), m.vert[c].cP());
        if (l > bestLength) { bestLength = l; best = c; }
      }
      return best;
    };

    std::vector<int> order(n, 0);
    std::vector<int> stemStart(n, -1); // first vertex of the stem each vertex lies on
    const int r = t.bfs.front();
    order[r] = 1;
    stemStart[r] = r;
    for (int v : t.bfs)
    {
      if (children[v].empty()) continue;

      int main = children[v].front();
      if (children[v].size() > 1)
      {
        const CoordType stemDir = m.vert[v].cP() - m.vert[stemStart[v]].cP();
        if (rule == LongestPath || stemDir.SquaredNorm() == 0)
          main = longestChild(v);
        else
        {
          ScalarType bestAngle = std::numeric_limits<ScalarType>::max();
          for (int c : children[v])
          {
            // vcg::Angle returns -1 for a zero-length vector, which would make a
            // degenerate child win; such a child cannot define a direction.
            const ScalarType a = Angle(stemDir, m.vert[c].cP() - m.vert[v].cP());
            if (a >= 0 && a < bestAngle) { bestAngle = a; main = c; }
          }
        }
      }

      for (int c : children[v])
      {
        if (c == main) { order[c] = order[v];     stemStart[c] = stemStart[v]; }
        else           { order[c] = order[v] + 1; stemStart[c] = v; }
      }
    }
    return order;
  }
};

} // end namespace tri
} // end namespace vcg

#endif // VCG_STREAM_ORDER_H
