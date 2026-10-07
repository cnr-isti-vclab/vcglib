/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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

#ifndef __VCGLIB_PLANAR_POLYGON_TESSELLATOR
#define __VCGLIB_PLANAR_POLYGON_TESSELLATOR

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <queue>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
namespace planar_polygon_detail {

/// The exact part of Orient2D, apart so that the filtered test stays small enough to be
/// inlined in the tessellators' inner loops; nearly collinear points are rare.
inline long double Orient2DExact(const Point2d &a, const Point2d &b, const Point2d &c)
{
	// Exact: each difference is hi + lo, each product of two parts is p + e, and the
	// sixteen terms are summed into a nonoverlapping expansion (increasing magnitude),
	// whose largest nonzero component carries the sign of the whole.
	auto twoDiff = [](double x, double y, double &lo) {
		const double s = x - y, bv = x - s;
		lo = (x - (s + bv)) + (bv - y);
		return s;
	};
	double dxb[2], dyc[2], dyb[2], dxc[2];
	dxb[0] = twoDiff(b.X(), a.X(), dxb[1]);
	dyc[0] = twoDiff(c.Y(), a.Y(), dyc[1]);
	dyb[0] = twoDiff(b.Y(), a.Y(), dyb[1]);
	dxc[0] = twoDiff(c.X(), a.X(), dxc[1]);
	double expansion[16];   // sixteen terms never need more components than that
	int size = 0;
	auto add = [&expansion, &size](double q) {   // Shewchuk's Grow-Expansion
		int n = 0;
		for (int i = 0; i < size; ++i) {
			const double h = expansion[i], s = q + h, bv = s - q;
			const double err = (q - (s - bv)) + (h - bv);
			if (err != 0) expansion[n++] = err;
			q = s;
		}
		size = n;
		if (q != 0) expansion[size++] = q;
	};
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j) {
			const double p = dxb[i] * dyc[j], m = dyb[i] * dxc[j];
			add(p); add(std::fma(dxb[i], dyc[j], -p));
			add(-m); add(-std::fma(dyb[i], dxc[j], -m));
		}
	return size == 0 ? 0.0L : static_cast<long double>(expansion[size - 1]);
}

/// Twice the signed area of the triangle (a, b, c): positive when it turns
/// counterclockwise, negative when clockwise, zero when the points are collinear.
///
/// The sign is exact (Shewchuk's adaptive orientation test): the determinant is
/// computed in double together with a bound on its rounding error, and when it does
/// not clear that bound -- nearly collinear points -- it is recomputed exactly, from
/// error-free differences and products (std::fma) summed as a floating-point expansion.
/// So callers comparing with 0 get decisions that agree with each other, which an
/// epsilon cannot promise; the value is only approximate. Exactness assumes IEEE
/// double arithmetic (no -ffast-math) and coordinates far from overflow and underflow.
/// (Before this, the plain determinant was computed in long double, which on many
/// platforms, Apple Silicon among them, is just double.)
inline long double Orient2D(const Point2d &a, const Point2d &b, const Point2d &c)
{
	const double l = (b.X() - a.X()) * (c.Y() - a.Y());
	const double r = (b.Y() - a.Y()) * (c.X() - a.X());
	const double det = l - r;
	const double eps = std::numeric_limits<double>::epsilon() / 2;   // unit roundoff
	if (std::abs(det) >= (3 + 16 * eps) * eps * (std::abs(l) + std::abs(r)))
		return det;
	return Orient2DExact(a, b, c);
}

inline bool PointInTriangle(
	const Point2d &point,
	const Point2d &a,
	const Point2d &b,
	const Point2d &c,
	long double winding,
	long double epsilon)
{
	return winding * Orient2D(a, b, point) >= -epsilon
		&& winding * Orient2D(b, c, point) >= -epsilon
		&& winding * Orient2D(c, a, point) >= -epsilon;
}

struct IndexedPoint2
{
	Point2d point;
	int index;
};

inline long double SignedDoubleArea(const std::vector<IndexedPoint2> &contour)
{
	long double area = 0;
	for (size_t i = 0; i < contour.size(); ++i) {
		const Point2d &a = contour[i].point;
		const Point2d &b = contour[(i + 1) % contour.size()].point;
		area += static_cast<long double>(a.X()) * static_cast<long double>(b.Y()) -
		        static_cast<long double>(a.Y()) * static_cast<long double>(b.X());
	}
	return area;
}

inline bool SamePoint(const Point2d &a, const Point2d &b, long double epsilon)
{
	const long double dx = static_cast<long double>(a.X()) - b.X();
	const long double dy = static_cast<long double>(a.Y()) - b.Y();
	return dx * dx + dy * dy <= epsilon * epsilon;
}

inline bool PointOnSegment(
	const Point2d &point,
	const Point2d &a,
	const Point2d &b,
	long double epsilon)
{
	if (std::abs(Orient2D(a, b, point)) > epsilon)
		return false;
	const long double edgeScale = std::max(
		std::abs(static_cast<long double>(b.X()) - a.X()),
		std::abs(static_cast<long double>(b.Y()) - a.Y()));
	const long double coordinateEpsilon = epsilon
		/ std::max(edgeScale, std::sqrt(epsilon));
	return point.X() >= std::min(a.X(), b.X()) - coordinateEpsilon
		&& point.X() <= std::max(a.X(), b.X()) + coordinateEpsilon
		&& point.Y() >= std::min(a.Y(), b.Y()) - coordinateEpsilon
		&& point.Y() <= std::max(a.Y(), b.Y()) + coordinateEpsilon;
}

inline int OrientationSign(long double value, long double epsilon)
{
	return value > epsilon ? 1 : (value < -epsilon ? -1 : 0);
}

inline bool SegmentsIntersect(
	const Point2d &a,
	const Point2d &b,
	const Point2d &c,
	const Point2d &d,
	long double epsilon)
{
	const int abc = OrientationSign(Orient2D(a, b, c), epsilon);
	const int abd = OrientationSign(Orient2D(a, b, d), epsilon);
	const int cda = OrientationSign(Orient2D(c, d, a), epsilon);
	const int cdb = OrientationSign(Orient2D(c, d, b), epsilon);
	if (abc * abd < 0 && cda * cdb < 0)
		return true;
	return (abc == 0 && PointOnSegment(c, a, b, epsilon))
		|| (abd == 0 && PointOnSegment(d, a, b, epsilon))
		|| (cda == 0 && PointOnSegment(a, c, d, epsilon))
		|| (cdb == 0 && PointOnSegment(b, c, d, epsilon));
}

inline bool PointInContour(
	const Point2d &point,
	const std::vector<IndexedPoint2> &contour,
	long double epsilon)
{
	bool inside = false;
	for (size_t i = 0, j = contour.size() - 1; i < contour.size(); j = i++) {
		const Point2d &a = contour[j].point;
		const Point2d &b = contour[i].point;
		if (PointOnSegment(point, a, b, epsilon))
			return false;
		if ((a.Y() > point.Y()) != (b.Y() > point.Y())) {
			const long double intersectionX = static_cast<long double>(a.X()) +
				(static_cast<long double>(b.X()) - a.X()) *
				(static_cast<long double>(point.Y()) - a.Y()) / (b.Y() - a.Y());
			if (static_cast<long double>(point.X()) < intersectionX)
				inside = !inside;
		}
	}
	return inside;
}

inline bool IsSimpleContour(
	const std::vector<IndexedPoint2> &contour,
	long double epsilon,
	long double pointEpsilon)
{
	const size_t count = contour.size();
	for (size_t i = 0; i < count; ++i) {
		if (SamePoint(contour[i].point, contour[(i + 1) % count].point, pointEpsilon))
			return false;
		for (size_t j = i + 1; j < count; ++j) {
			if (SamePoint(contour[i].point, contour[j].point, pointEpsilon))
				return false;
		}
	}
	for (size_t i = 0; i < count; ++i) {
		const size_t iNext = (i + 1) % count;
		for (size_t j = i + 1; j < count; ++j) {
			const size_t jNext = (j + 1) % count;
			if (i == j || iNext == j || jNext == i)
				continue;
			if (SegmentsIntersect(
					contour[i].point, contour[iNext].point,
					contour[j].point, contour[jNext].point, epsilon))
				return false;
		}
	}
	return true;
}

inline bool ContoursIntersect(
	const std::vector<IndexedPoint2> &first,
	const std::vector<IndexedPoint2> &second,
	long double epsilon)
{
	for (size_t i = 0; i < first.size(); ++i)
		for (size_t j = 0; j < second.size(); ++j)
			if (SegmentsIntersect(
					first[i].point, first[(i + 1) % first.size()].point,
					second[j].point, second[(j + 1) % second.size()].point,
					epsilon))
				return true;
	return false;
}

inline bool BridgeIsVisible(
	const std::vector<IndexedPoint2> &merged,
	size_t outerIndex,
	const std::vector<IndexedPoint2> &hole,
	size_t holeIndex,
	const std::vector<std::vector<IndexedPoint2>> &allHoles,
	long double epsilon)
{
	const Point2d &a = merged[outerIndex].point;
	const Point2d &b = hole[holeIndex].point;
	for (size_t i = 0; i < merged.size(); ++i) {
		const size_t next = (i + 1) % merged.size();
		if (merged[i].index == merged[outerIndex].index
			|| merged[next].index == merged[outerIndex].index)
			continue;
		if (SegmentsIntersect(a, b, merged[i].point, merged[next].point, epsilon))
			return false;
	}
	for (const auto &candidateHole : allHoles) {
		for (size_t i = 0; i < candidateHole.size(); ++i) {
			const size_t next = (i + 1) % candidateHole.size();
			const bool isHoleEndpointEdge = &candidateHole == &hole
				&& (i == holeIndex || next == holeIndex);
			const bool isMergedEndpointEdge =
				candidateHole[i].index == merged[outerIndex].index
				|| candidateHole[next].index == merged[outerIndex].index;
			if (!isHoleEndpointEdge && !isMergedEndpointEdge && SegmentsIntersect(
					a, b, candidateHole[i].point, candidateHole[next].point, epsilon))
				return false;
		}
	}
	const Point2d midpoint((a.X() + b.X()) * 0.5, (a.Y() + b.Y()) * 0.5);
	if (!PointInContour(midpoint, merged, epsilon))
		return false;
	for (const auto &candidateHole : allHoles)
		if (PointInContour(midpoint, candidateHole, epsilon))
			return false;
	return true;
}

inline bool MergeHole(
	std::vector<IndexedPoint2> &merged,
	const std::vector<IndexedPoint2> &hole,
	const std::vector<std::vector<IndexedPoint2>> &allHoles,
	long double epsilon)
{
	size_t bestOuter = 0;
	size_t bestHole = 0;
	long double bestLength2 = std::numeric_limits<long double>::max();
	bool found = false;
	for (size_t i = 0; i < merged.size(); ++i) {
		for (size_t j = 0; j < hole.size(); ++j) {
			if (!BridgeIsVisible(merged, i, hole, j, allHoles, epsilon))
				continue;
			const long double dx = static_cast<long double>(merged[i].point.X()) - hole[j].point.X();
			const long double dy = static_cast<long double>(merged[i].point.Y()) - hole[j].point.Y();
			const long double length2 = dx * dx + dy * dy;
			if (length2 < bestLength2) {
				bestLength2 = length2;
				bestOuter = i;
				bestHole = j;
				found = true;
			}
		}
	}
	if (!found)
		return false;

	std::vector<IndexedPoint2> result;
	result.reserve(merged.size() + hole.size() + 2);
	result.insert(result.end(), merged.begin(), merged.begin() + bestOuter + 1);
	for (size_t i = 0; i < hole.size(); ++i)
		result.push_back(hole[(bestHole + i) % hole.size()]);
	result.push_back(hole[bestHole]);
	result.push_back(merged[bestOuter]);
	result.insert(result.end(), merged.begin() + bestOuter + 1, merged.end());
	merged.swap(result);
	return true;
}

inline bool TessellateWeaklySimpleContour(
	const std::vector<IndexedPoint2> &contour,
	std::vector<int> &triangles,
	long double epsilon)
{
	const size_t count = contour.size();
	std::vector<int> previous(count), next(count);
	std::vector<unsigned char> active(count, 1);
	for (size_t i = 0; i < count; ++i) {
		previous[i] = int((i + count - 1) % count);
		next[i] = int((i + 1) % count);
	}

	int current = 0;
	int activeCount = int(count);
	while (activeCount > 2) {
		bool foundEar = false;
		for (int attempts = 0; attempts < activeCount; ++attempts) {
			const int before = previous[current];
			const int after = next[current];
			if (Orient2D(contour[before].point, contour[current].point, contour[after].point) > epsilon) {
				bool containsVertex = false;
				for (size_t candidate = 0; candidate < count; ++candidate) {
					if (!active[candidate] || int(candidate) == before
						|| int(candidate) == current || int(candidate) == after)
						continue;
					if (contour[candidate].index == contour[before].index
						|| contour[candidate].index == contour[current].index
						|| contour[candidate].index == contour[after].index)
						continue;
					if (PointInTriangle(
							contour[candidate].point,
							contour[before].point,
							contour[current].point,
							contour[after].point,
							1, epsilon)) {
						containsVertex = true;
						break;
					}
				}
				if (!containsVertex) {
					triangles.push_back(contour[before].index);
					triangles.push_back(contour[current].index);
					triangles.push_back(contour[after].index);
					next[before] = after;
					previous[after] = before;
					active[current] = 0;
					current = after;
					--activeCount;
					foundEar = true;
					break;
				}
			}
			current = next[current];
		}
		if (!foundEar)
			return false;
	}
	return true;
}

} // namespace planar_polygon_detail

/**
 * Triangulate one finite, simple 2D polygon with ear clipping.
 *
 * The input must be one boundary loop, without holes, repeated consecutive
 * points, or self-intersections. Output indices are local to points, retain
 * its winding, and are appended only on success. Collinear vertices are supported;
 * no zero-area triangles are emitted. False means the polygon is degenerate,
 * invalid, or could not be triangulated reliably, and leaves output unchanged.
 */
template <class POINT_CONTAINER>
bool TessellatePlanarPolygon2(const POINT_CONTAINER &points, std::vector<int> &output)
{
	using namespace planar_polygon_detail;
	const size_t count = points.size();
	if (count < 3)
		return false;

	// Work in translated double coordinates. Translation avoids losing the small
	// differences that define the polygon when coordinates have a large offset.
	std::vector<Point2d> p(count);
	const double originX = double(points[0][0]);
	const double originY = double(points[0][1]);
	double scale = 0;
	for (size_t i = 0; i < count; ++i) {
		const double x = double(points[i][0]) - originX;
		const double y = double(points[i][1]) - originY;
		if (!std::isfinite(x) || !std::isfinite(y))
			return false;
		p[i] = Point2d(x, y);
		scale = std::max(scale, std::max(std::abs(x), std::abs(y)));
	}
	if (!(scale > 0))
		return false;

	const long double scale2 = static_cast<long double>(scale) * scale;
	const long double epsilon = scale2 * 64 * std::numeric_limits<double>::epsilon();
	const long double edgeEpsilon2 = scale2 * 1024
		* std::numeric_limits<double>::epsilon() * std::numeric_limits<double>::epsilon();

	long double signedDoubleArea = 0;
	for (size_t i = 0; i < count; ++i) {
		const Point2d &a = p[i];
		const Point2d &b = p[(i + 1) % count];
		const long double dx = static_cast<long double>(b.X()) - static_cast<long double>(a.X());
		const long double dy = static_cast<long double>(b.Y()) - static_cast<long double>(a.Y());
		if (dx * dx + dy * dy <= edgeEpsilon2)
			return false;
		signedDoubleArea += static_cast<long double>(a.X()) * static_cast<long double>(b.Y()) -
		                    static_cast<long double>(a.Y()) * static_cast<long double>(b.X());
	}
	if (std::abs(signedDoubleArea) <= epsilon * count)
		return false;
	const long double winding = signedDoubleArea > 0 ? 1 : -1;

	std::vector<int> previous(count), next(count);
	std::vector<unsigned char> active(count, 1);
	for (size_t i = 0; i < count; ++i) {
		previous[i] = int((i + count - 1) % count);
		next[i] = int((i + 1) % count);
	}

	std::vector<int> triangles;
	triangles.reserve(3 * (count - 2));
	int current = 0;
	int activeCount = int(count);
	while (activeCount > 2) {
		bool foundEar = false;
		for (int attempts = 0; attempts < activeCount; ++attempts) {
			const int before = previous[current];
			const int after = next[current];
			if (winding * Orient2D(p[before], p[current], p[after]) > epsilon) {
				// Convexity alone is insufficient: another active vertex inside
				// this triangle means the replacement diagonal is not an ear.
				bool containsVertex = false;
				for (size_t candidate = 0; candidate < count; ++candidate) {
					if (!active[candidate] || int(candidate) == before
						|| int(candidate) == current || int(candidate) == after)
						continue;
					if (PointInTriangle(
							p[candidate], p[before], p[current], p[after], winding, epsilon)) {
						containsVertex = true;
						break;
					}
				}
				if (!containsVertex) {
					triangles.push_back(before);
					triangles.push_back(current);
					triangles.push_back(after);
					next[before] = after;
					previous[after] = before;
					active[current] = 0;
					current = after;
					--activeCount;
					foundEar = true;
					break;
				}
			}
			current = next[current];
		}
		if (!foundEar)
			return false;
	}

	if (triangles.size() != 3 * (count - 2))
		return false;

	// Validate the postcondition independently from ear selection. Consistent
	// winding and equal area reject overlaps, outside triangles, and gaps.
	long double triangleDoubleArea = 0;
	for (size_t i = 0; i < triangles.size(); i += 3) {
		const long double area = winding * Orient2D(
			p[triangles[i]], p[triangles[i + 1]], p[triangles[i + 2]]);
		if (area <= epsilon)
			return false;
		triangleDoubleArea += area;
	}
	const long double areaTolerance = std::max(
		std::abs(signedDoubleArea) * 1e-12L, epsilon * count * 4);
	if (std::abs(triangleDoubleArea - std::abs(signedDoubleArea)) > areaTolerance)
		return false;

	output.insert(output.end(), triangles.begin(), triangles.end());
	return true;
}

namespace planar_polygon_detail {

/// A triangulation as counterclockwise index triples, with the map from each directed
/// edge to its triangle that flips and insertions need. Shared by the constrained Delaunay
/// flips and the quality refinement. \a orient decides convexity and should have an exact
/// sign; \a len2 is the squared length angles are measured with, which may differ from the
/// frame \a orient works in (a projection or a barycentric frame keeps orientations, not
/// angles); \a fixed(a, b) says whether an edge must stay.
///
/// Optionally journaled: between BeginJournal() and Undo() every triangle written is
/// remembered, so an insertion that turns out to be unwanted can be taken back.
template <class ORIENT, class LEN2, class FIXED>
struct FlipTriangulation
{
	std::vector<int> &tris;
	ORIENT orient;
	LEN2 len2;
	FIXED fixed;
	std::unordered_map<std::uint64_t, size_t> corner;   // directed edge (a, b) -> position of a in tris

	FlipTriangulation(std::vector<int> &t, ORIENT o, LEN2 l, FIXED f) : tris(t), orient(o), len2(l), fixed(f)
	{
		for (size_t i = 0; i < tris.size(); i += 3) Index(i);
	}
	static std::uint64_t Key(int a, int b) { return (std::uint64_t(std::uint32_t(a)) << 32) | std::uint32_t(b); }
	void Index(size_t t) { for (int k = 0; k < 3; ++k) corner[Key(tris[t + k], tris[t + (k + 1) % 3])] = t + k; }
	/// Forget the edges of the triangle at t -- those still mapped to it: rewriting two
	/// triangles in turn, as a flip does, hands an edge from one to the other in between.
	void Unindex(size_t t)
	{
		for (int k = 0; k < 3; ++k) {
			const auto it = corner.find(Key(tris[t + k], tris[t + (k + 1) % 3]));
			if (it != corner.end() && it->second - it->second % 3 == t) corner.erase(it);
		}
	}
	/// The position in tris of a, in the triangle that has the directed edge (a, b).
	bool Find(int a, int b, size_t &pos) const
	{
		const auto it = corner.find(Key(a, b));
		if (it == corner.end()) return false;
		pos = it->second;
		return true;
	}
	int Next(size_t pos) const { return tris[pos - pos % 3 + (pos + 1) % 3]; }
	int Apex(size_t pos) const { return tris[pos - pos % 3 + (pos + 2) % 3]; }

	// The journal: the triangles there were before, and how many.
	bool journaling = false;
	size_t journalSize = 0;
	std::vector<std::pair<size_t, std::array<int, 3>>> journal;
	void BeginJournal() { journaling = true; journalSize = tris.size(); journal.clear(); }
	void EndJournal() { journaling = false; }
	void Set(size_t t, int a, int b, int c)
	{
		if (journaling && t < journalSize
		    && std::find_if(journal.begin(), journal.end(), [t](const auto &j) { return j.first == t; }) == journal.end())
			journal.push_back({t, {{tris[t], tris[t + 1], tris[t + 2]}}});
		Unindex(t);
		tris[t] = a; tris[t + 1] = b; tris[t + 2] = c;
		Index(t);
	}
	size_t Add(int a, int b, int c)
	{
		tris.insert(tris.end(), { a, b, c });
		Index(tris.size() - 3);
		return tris.size() - 3;
	}
	void Undo()
	{
		for (size_t t = journalSize; t < tris.size(); t += 3) Unindex(t);
		tris.resize(journalSize);
		for (const auto &j : journal) { Unindex(j.first); for (int k = 0; k < 3; ++k) tris[j.first + k] = j.second[k]; }
		for (const auto &j : journal) Index(j.first);
		journaling = false;
	}

	/// The largest cosine among the angles of (a, b, c): larger as its smallest angle shrinks.
	double MaxCos(int a, int b, int c) const
	{
		const double ab = len2(a, b), bc = len2(b, c), ca = len2(c, a);
		return std::max({ (ab + ca - bc) / (2 * std::sqrt(ab * ca)),
		                  (ab + bc - ca) / (2 * std::sqrt(ab * bc)),
		                  (bc + ca - ab) / (2 * std::sqrt(bc * ca)) });
	}
	/// Lawson's flips over the queued edges: an edge that is not fixed is flipped when its two
	/// triangles form a strictly convex quad and the flip increases the smaller of their
	/// smallest angles, which for a convex quad is the Delaunay choice. Each flip strictly
	/// improves the sorted angles, so the flips stop. Boundary edges, with one triangle, stay.
	void Legalize(std::deque<std::pair<int, int>> &queue)
	{
		while (!queue.empty()) {
			const int a = queue.front().first, b = queue.front().second;
			queue.pop_front();
			size_t p1, p2;
			if (fixed(a, b) || !Find(a, b, p1) || !Find(b, a, p2)) continue;
			const int c = Apex(p1), d = Apex(p2);
			if (!(orient(c, a, d) > 0 && orient(d, b, c) > 0)) continue;
			if (!(std::max(MaxCos(c, a, d), MaxCos(d, b, c)) < std::max(MaxCos(a, b, c), MaxCos(b, a, d)) - 1e-12)) continue;
			Set(p1 - p1 % 3, c, a, d);
			Set(p2 - p2 % 3, d, b, c);
			queue.insert(queue.end(), { {a, d}, {d, b}, {b, c}, {c, a} });
		}
	}
	void LegalizeAll()
	{
		std::deque<std::pair<int, int>> queue;
		for (size_t i = 0; i < tris.size(); ++i)
			if (tris[i] < Next(i)) queue.push_back({tris[i], Next(i)});
		Legalize(queue);
	}
};

template <class ORIENT, class LEN2, class FIXED>
FlipTriangulation<ORIENT, LEN2, FIXED> MakeFlipTriangulation(std::vector<int> &t, ORIENT o, LEN2 l, FIXED f)
{
	return FlipTriangulation<ORIENT, LEN2, FIXED>(t, o, l, f);
}

} // namespace planar_polygon_detail

/**
 * Flip a triangulation to the constrained Delaunay one (Lawson's flips).
 *
 * \a tris holds counterclockwise triangles as index triples. An interior edge that is
 * not in \a fixed (pairs, smaller index first) is flipped when its two triangles form a
 * strictly convex quad and the flip increases the smaller of their smallest angles, which
 * for a convex quad is the Delaunay choice; each flip strictly improves the sorted angles,
 * so the flips stop. Boundary edges, having one triangle, never flip. \a orient(a, b, c)
 * decides convexity and should have an exact sign (Orient2D); \a len2(a, b) is the squared
 * length the angles are measured with, which may differ from the frame \a orient works in:
 * a projection or a barycentric frame keeps orientations but not angles.
 *
 * Used where triangle quality matters: the faces rebuilt around an embedded curve
 * (CoMEmbed) and, on request, the planar contour tessellators.
 */
template <class ORIENT, class LEN2>
void FlipToConstrainedDelaunay(std::vector<int> &tris, const std::set<std::pair<int, int>> &fixed,
                               ORIENT orient, LEN2 len2)
{
	auto ft = planar_polygon_detail::MakeFlipTriangulation(tris, orient, len2,
		[&fixed](int a, int b) { return fixed.count({std::min(a, b), std::max(a, b)}) > 0; });
	ft.LegalizeAll();
}

/// How RefinePlanarTriangulation2 refines.
struct PlanarRefinement
{
	double minAngle = 20;        ///< degrees; no triangle should have a smaller angle, unless the input forces it
	bool splitBoundary = false;  ///< may points be added on the boundary edges? Needed for the guarantee near the boundary
	int maxSteiner = -1;         ///< at most this many points; negative: twenty times the input points, plus a hundred
};

/// A point RefinePlanarTriangulation2 added, and how: on the boundary edge from input point
/// \a a to input point \a b at parameter \a t (\a c < 0), or inside the triangle (\a a, \a b,
/// \a c) -- input or earlier added points -- with barycentric weights \a w. Enough for a
/// caller to give it the attributes of where it lies, in the order the points were added.
struct SteinerPoint
{
	int a = -1, b = -1, c = -1;
	double t = 0;
	double w[3] = { 0, 0, 0 };
};

/**
 * Refine a constrained Delaunay triangulation of a planar region by adding points, until
 * no triangle has an angle below \a opt.minAngle (Delaunay refinement: Ruppert 1995, with
 * the improvements described by Shewchuk 2002).
 *
 * \a points are the vertices, \a tris counterclockwise triples, Delaunay with respect to
 * their boundary edges (FlipToConstrainedDelaunay, or TessellatePlanarContours2 with
 * delaunay set). The boundary is the edges that have one triangle; it is the only
 * constraint. Added points are appended to \a points, the triangulation is rewritten in
 * \a tris, and each added point is described in \a added.
 *
 * A bad triangle is split at its circumcenter, or at Ungor's off-center, closer in, when
 * that is enough: the apex of the isosceles triangle on its shortest edge whose apex angle
 * is the target, which reaches the same quality with fewer points. A point that would
 * encroach on a boundary edge (lie inside the circle having that edge as diameter) is not
 * inserted; the edge is split instead, at its midpoint or, next to an input vertex, at a
 * power-of-two distance from it ("concentric shells"), so two edges meeting at a small
 * angle do not split each other forever. A thin triangle whose small angle comes from the
 * input itself -- its shortest edge joins two points on two boundary edges, equidistant
 * from the input vertex where those meet -- is left alone (Miller, Pav and Walkington).
 *
 * Without \a opt.splitBoundary the boundary is kept as it is: a point that would land
 * beyond a boundary edge, or encroach on one, is not inserted, and the triangles along the
 * boundary are improved only as far as that allows. Points on a boundary edge are on it exactly by construction: tests
 * involving them and that edge treat them as collinear.
 *
 * Returns the number of points added. Below about 20.7 degrees the refinement provably
 * stops; up to about 33 it does in practice; the budget \a opt.maxSteiner stops it in any
 * case, and \a opt.minAngle is clamped to 34.
 */
inline int RefinePlanarTriangulation2(std::vector<Point2d> &points, std::vector<int> &tris,
                                      const PlanarRefinement &opt, std::vector<SteinerPoint> &added)
{
	using namespace planar_polygon_detail;
	const int inputCount = int(points.size());
	const double pi = std::acos(-1.0);
	const double minAngle = std::min(opt.minAngle, 34.0) * pi / 180;
	if (!(minAngle > 0) || tris.empty()) return 0;
	const double cosMax = std::cos(minAngle);
	const double offDistance = 0.475 / std::tan(minAngle / 2);   // from the shortest edge, times its length
	const int budget = opt.maxSteiner >= 0 ? opt.maxSteiner : 20 * inputCount + 100;
	typedef std::uint64_t Key;
	const auto ukey = [](int a, int b) { return (Key(std::uint32_t(std::min(a, b))) << 32) | std::uint32_t(std::max(a, b)); };

	// The boundary edges, each remembered as a piece of the input edge it came from.
	std::vector<std::array<int, 2>> segments;            // the input boundary edges
	std::unordered_map<Key, int> segmentOf;               // current boundary edge -> input edge
	std::vector<int> pointSegment(points.size(), -1);     // the input edge a point was added on
	{
		std::unordered_map<Key, int> count;
		for (size_t i = 0; i < tris.size(); ++i) ++count[ukey(tris[i], tris[i - i % 3 + (i + 1) % 3])];
		for (size_t i = 0; i < tris.size(); ++i) {
			const int a = tris[i], b = tris[i - i % 3 + (i + 1) % 3];
			if (count[ukey(a, b)] == 1) { segmentOf[ukey(a, b)] = int(segments.size()); segments.push_back({{a, b}}); }
		}
	}
	const auto onSegment = [&](int p, int s) {
		return p == segments[size_t(s)][0] || p == segments[size_t(s)][1] || pointSegment[size_t(p)] == s;
	};
	// Points added on an input edge are on it exactly: three points of one edge are collinear.
	const auto orient = [&](int a, int b, int c) -> long double {
		for (int s : { pointSegment[size_t(a)], pointSegment[size_t(b)], pointSegment[size_t(c)] })
			if (s >= 0 && onSegment(a, s) && onSegment(b, s) && onSegment(c, s)) return 0;
		return Orient2D(points[size_t(a)], points[size_t(b)], points[size_t(c)]);
	};
	const auto len2 = [&](int a, int b) { return (points[size_t(a)] - points[size_t(b)]).SquaredNorm(); };
	const auto isBoundary = [&](int a, int b) { return segmentOf.count(ukey(a, b)) > 0; };
	auto ft = MakeFlipTriangulation(tris, orient, len2, isBoundary);
	// p inside the circle that has ab as diameter.
	const auto encroaches = [&](int p, int a, int b) {
		return (points[size_t(a)] - points[size_t(p)]).dot(points[size_t(b)] - points[size_t(p)]) < 0;
	};

	// The bad triangles, worst first: keyed by the cosine of the smallest angle.
	std::priority_queue<std::pair<double, std::array<int, 3>>> bad;
	const auto test = [&](int a, int b, int c) {
		if (!(len2(a, b) > 0 && len2(b, c) > 0 && len2(c, a) > 0)) return;
		const double cs = ft.MaxCos(a, b, c);
		if (!(cs > cosMax)) return;
		int p = a, q = b;   // the shortest edge
		if (len2(b, c) < len2(p, q)) { p = b; q = c; }
		if (len2(c, a) < len2(p, q)) { p = c; q = a; }
		// The small angle is the input's: the shortest edge joins points added on two input
		// edges, equidistant from the input vertex where these meet. Splitting would not end.
		const int sp = pointSegment[size_t(p)], sq = pointSegment[size_t(q)];
		if (sp >= 0 && sq >= 0 && sp != sq)
			for (int j : segments[size_t(sp)])
				if (j == segments[size_t(sq)][0] || j == segments[size_t(sq)][1]) {
					const double d1 = len2(j, p), d2 = len2(j, q);
					if (d1 < 1.001 * d2 && d1 > 0.999 * d2) return;
				}
		bad.push({cs, {{a, b, c}}});
	};
	for (size_t t = 0; t < tris.size(); t += 3) test(tris[t], tris[t + 1], tris[t + 2]);

	// The corners of v in the triangles around it, starting from the triangle at slot t.
	const auto star = [&](int v, size_t t) {
		size_t start = t;
		while (tris[start] != v) ++start;
		std::vector<size_t> around;
		size_t cur = start, next;
		for (;;) {   // counterclockwise: the next triangle has the directed edge (v, apex)
			around.push_back(cur);
			if (!ft.Find(v, ft.Apex(cur), next)) break;   // the boundary
			if (next == start) return around;            // all the way round
			cur = next;
		}
		for (cur = start;;) {   // stopped by the boundary: clockwise from the start as well
			size_t back;
			if (!ft.Find(ft.Next(cur), v, back)) break;
			cur = back - back % 3 + (back % 3 + 1) % 3;   // v's corner in that triangle
			around.push_back(cur);
		}
		return around;
	};
	std::deque<std::pair<int, int>> encroached;   // boundary edges to split
	// After v went in: test the triangles around it, and list the boundary edges it encroaches.
	const auto afterInsert = [&](int v, size_t t, std::vector<std::pair<int, int>> *enc) {
		for (size_t c : star(v, t)) {
			const int n = ft.Next(c), o = ft.Apex(c);
			if (isBoundary(n, o) && encroaches(v, n, o)) {
				if (enc) enc->push_back({n, o});
			}
			test(v, n, o);
		}
	};
	if (opt.splitBoundary)
		for (size_t i = 0; i < tris.size(); ++i)
			if (isBoundary(tris[i], ft.Next(i)) && encroaches(ft.Apex(i), tris[i], ft.Next(i)))
				encroached.push_back({tris[i], ft.Next(i)});

	int steiner = 0;
	// Split the boundary edge (a, b), triangle (a, b, c), at the point x on it.
	const auto splitBoundaryEdge = [&](int a, int b) {
		size_t pos;
		if (!ft.Find(a, b, pos)) return;
		const int s = segmentOf[ukey(a, b)];
		const int c = ft.Apex(pos);
		// Midway, or next to an input vertex at a power-of-two distance from it: edges that
		// meet at a small angle are then split at the same radii and stop encroaching.
		const double l = std::sqrt(len2(a, b));
		double f = 0.5;
		if ((a < inputCount) != (b < inputCount)) {
			double p2 = 1;
			while (l > 3 * p2) p2 *= 2;
			while (l < 1.5 * p2) p2 /= 2;
			f = (a < inputCount) ? p2 / l : 1 - p2 / l;
		}
		const int x = int(points.size());
		points.push_back(points[size_t(a)] + (points[size_t(b)] - points[size_t(a)]) * f);
		pointSegment.push_back(s);
		segmentOf.erase(ukey(a, b));
		segmentOf[ukey(a, x)] = s;
		segmentOf[ukey(x, b)] = s;
		const size_t t = pos - pos % 3;
		ft.Set(t, a, x, c);
		ft.Add(x, b, c);
		std::deque<std::pair<int, int>> q{ {b, c}, {c, a} };
		ft.Legalize(q);
		SteinerPoint sp;
		const int o0 = segments[size_t(s)][0], o1 = segments[size_t(s)][1];
		sp.a = o0; sp.b = o1;
		sp.t = std::sqrt(len2(o0, x) / len2(o0, o1));
		added.push_back(sp);
		++steiner;
		std::vector<std::pair<int, int>> enc;
		afterInsert(x, t, &enc);
		encroached.insert(encroached.end(), enc.begin(), enc.end());
	};

	while (steiner < budget) {
		if (!encroached.empty()) {
			// Split whatever is still a boundary edge: it was encroached on by a vertex, by a
			// point that was taken back, or lay between a bad triangle and its circumcenter.
			const auto e = encroached.front();
			encroached.pop_front();
			if (isBoundary(e.first, e.second)) splitBoundaryEdge(e.first, e.second);
			continue;
		}
		if (bad.empty()) break;
		const std::array<int, 3> tri = bad.top().second;
		bad.pop();
		size_t pos;
		if (!ft.Find(tri[0], tri[1], pos) || ft.Apex(pos) != tri[2]) continue;   // gone meanwhile
		// The circumcenter, or the off-center when it is closer to the shortest edge.
		const Point2d &A = points[size_t(tri[0])], &B = points[size_t(tri[1])], &C = points[size_t(tri[2])];
		const Point2d ab = B - A, ac = C - A;
		const double d = 2 * (ab[0] * ac[1] - ab[1] * ac[0]);
		if (!(std::abs(d) > 0)) continue;
		const Point2d cc = A + Point2d((ac[1] * ab.SquaredNorm() - ab[1] * ac.SquaredNorm()) / d,
		                               (ab[0] * ac.SquaredNorm() - ac[0] * ab.SquaredNorm()) / d);
		int p = tri[0], q = tri[1];
		if (len2(tri[1], tri[2]) < len2(p, q)) { p = tri[1]; q = tri[2]; }
		if (len2(tri[2], tri[0]) < len2(p, q)) { p = tri[2]; q = tri[0]; }
		const Point2d mid = (points[size_t(p)] + points[size_t(q)]) * 0.5;
		const double h = offDistance * std::sqrt(len2(p, q)), toCc = (cc - mid).Norm();
		const Point2d x = (h < toCc) ? mid + (cc - mid) * (h / toCc) : cc;
		if (!std::isfinite(x[0]) || !std::isfinite(x[1])) continue;

		// Locate it, walking from the bad triangle.
		const int xi = int(points.size());
		points.push_back(x);
		pointSegment.push_back(-1);
		enum { Inside, OnEdge, OnBoundary, Beyond, Nowhere } where = Nowhere;
		size_t t = pos - pos % 3;
		int ea = -1, eb = -1;
		for (size_t step = 0; step < tris.size(); ++step) {
			int neg = -1, zeros = 0, zeroEdge = -1;
			for (int k = 0; k < 3; ++k) {
				const int k1 = int((k + step) % 3);
				const long double o = orient(tris[t + k1], tris[t + (k1 + 1) % 3], xi);
				if (o < 0 && neg < 0) neg = k1;
				if (o == 0) { ++zeros; zeroEdge = k1; }
			}
			if (neg >= 0) {
				ea = tris[t + neg]; eb = tris[t + (neg + 1) % 3];
				if (isBoundary(ea, eb)) { where = Beyond; break; }
				size_t back;
				if (!ft.Find(eb, ea, back)) break;
				t = back - back % 3;
				continue;
			}
			if (zeros == 0) where = Inside;
			else if (zeros == 1) {
				ea = tris[t + zeroEdge]; eb = tris[t + (zeroEdge + 1) % 3];
				where = isBoundary(ea, eb) ? OnBoundary : OnEdge;
			}
			break;
		}
		if (where != Inside && where != OnEdge) {
			points.pop_back();
			pointSegment.pop_back();
			// Beyond a boundary edge, or on one: that edge is in the way. Split it, if allowed,
			// and come back to this triangle.
			if (opt.splitBoundary && (where == Beyond || where == OnBoundary)) {
				encroached.push_back({ea, eb});
				bad.push({ft.MaxCos(tri[0], tri[1], tri[2]), tri});
			}
			continue;
		}
		// Where it goes, for the caller: the triangle it falls in and its weights there.
		SteinerPoint sp;
		sp.a = tris[t]; sp.b = tris[t + 1]; sp.c = tris[t + 2];
		{
			const long double area = Orient2D(points[size_t(sp.a)], points[size_t(sp.b)], points[size_t(sp.c)]);
			sp.w[0] = double(Orient2D(points[size_t(sp.b)], points[size_t(sp.c)], x) / area);
			sp.w[1] = double(Orient2D(points[size_t(sp.c)], points[size_t(sp.a)], x) / area);
			sp.w[2] = 1 - sp.w[0] - sp.w[1];
		}
		ft.BeginJournal();
		std::deque<std::pair<int, int>> legal;
		if (where == Inside) {
			const int a = tris[t], b = tris[t + 1], c = tris[t + 2];
			ft.Set(t, a, b, xi);
			ft.Add(b, c, xi);
			ft.Add(c, a, xi);
			legal = { {a, b}, {b, c}, {c, a} };
		} else {
			size_t p1, p2;
			ft.Find(ea, eb, p1);
			const bool twin = ft.Find(eb, ea, p2);
			const int c = ft.Apex(p1);
			const size_t t1 = p1 - p1 % 3;
			ft.Set(t1, ea, xi, c);
			ft.Add(xi, eb, c);
			legal = { {eb, c}, {c, ea} };
			if (twin) {
				const int dd = ft.Apex(p2);
				const size_t t2 = p2 - p2 % 3;
				ft.Set(t2, eb, xi, dd);
				ft.Add(xi, ea, dd);
				legal.insert(legal.end(), { {ea, dd}, {dd, eb} });
			}
			t = t1;
		}
		ft.Legalize(legal);
		std::vector<std::pair<int, int>> enc;
		for (size_t c : star(xi, t)) {
			const int n = ft.Next(c), o = ft.Apex(c);
			if (isBoundary(n, o) && encroaches(xi, n, o)) enc.push_back({n, o});
		}
		if (!enc.empty()) {
			// It would encroach on the boundary: take it back. The edges it encroaches on are
			// split instead, when that is allowed, and the triangle tried again; otherwise the
			// triangle stays as it is, since a point that close to a fixed edge would only
			// make thinner triangles against it.
			ft.Undo();
			points.pop_back();
			pointSegment.pop_back();
			if (opt.splitBoundary) {
				encroached.insert(encroached.end(), enc.begin(), enc.end());
				bad.push({ft.MaxCos(tri[0], tri[1], tri[2]), tri});
			}
			continue;
		}
		ft.EndJournal();
		added.push_back(sp);
		++steiner;
		afterInsert(xi, t, nullptr);
	}
	return steiner;
}

/**
 * Triangulate finite, simple 2D contours using the even-odd fill rule.
 *
 * Contours may have arbitrary winding and order and may describe disconnected
 * regions, holes, and islands nested to any depth. They must not cross, touch,
 * overlap, or contain repeated points. Output indices address the input points
 * flattened in contour order and all triangles are counter-clockwise. False
 * means that the input is invalid or cannot be triangulated reliably and leaves
 * output unchanged.
 *
 * Hole elimination follows the bridge-and-ear-clipping approach popularized by
 * FIST and Mapbox earcut.hpp. This implementation deliberately favors validation
 * and a small vcglib-style API over recovery of malformed polygon data.
 */
template <class CONTOUR_CONTAINER>
bool TessellatePlanarContours2(
	const CONTOUR_CONTAINER &inputContours,
	std::vector<int> &output,
	bool delaunay = false)
{
	using namespace planar_polygon_detail;
	if (inputContours.empty())
		return false;

	double originX = 0;
	double originY = 0;
	bool haveOrigin = false;
	double scale = 0;
	for (const auto &contour : inputContours) {
		if (contour.size() < 3)
			return false;
		for (const auto &point : contour) {
			const double x = double(point[0]);
			const double y = double(point[1]);
			if (!std::isfinite(x) || !std::isfinite(y))
				return false;
			if (!haveOrigin) {
				originX = x;
				originY = y;
				haveOrigin = true;
			}
			scale = std::max(scale, std::max(std::abs(x - originX), std::abs(y - originY)));
		}
	}
	if (!(scale > 0))
		return false;

	const long double scale2 = static_cast<long double>(scale) * scale;
	const long double epsilon = scale2 * 64 * std::numeric_limits<double>::epsilon();
	const long double pointEpsilon = static_cast<long double>(scale) * 32
		* std::numeric_limits<double>::epsilon();

	std::vector<std::vector<IndexedPoint2>> contours;
	std::vector<Point2d> flatPoints;
	contours.reserve(inputContours.size());
	int flatIndex = 0;
	for (const auto &inputContour : inputContours) {
		std::vector<IndexedPoint2> contour;
		contour.reserve(inputContour.size());
		for (const auto &point : inputContour) {
			const Point2d translated(double(point[0]) - originX, double(point[1]) - originY);
			contour.push_back({translated, flatIndex++});
			flatPoints.push_back(translated);
		}
		if (!IsSimpleContour(contour, epsilon, pointEpsilon))
			return false;
		contours.push_back(std::move(contour));
	}

	const size_t contourCount = contours.size();
	std::vector<long double> areas(contourCount);
	for (size_t i = 0; i < contourCount; ++i) {
		areas[i] = SignedDoubleArea(contours[i]);
		if (std::abs(areas[i]) <= epsilon * contours[i].size())
			return false;
		for (size_t j = 0; j < i; ++j)
			if (ContoursIntersect(contours[i], contours[j], epsilon))
				return false;
	}

	std::vector<int> parent(contourCount, -1);
	for (size_t i = 0; i < contourCount; ++i) {
		long double parentArea = std::numeric_limits<long double>::max();
		for (size_t j = 0; j < contourCount; ++j) {
			if (i == j || std::abs(areas[j]) <= std::abs(areas[i]))
				continue;
			if (PointInContour(contours[i][0].point, contours[j], epsilon)
				&& std::abs(areas[j]) < parentArea) {
				parent[i] = int(j);
				parentArea = std::abs(areas[j]);
			}
		}
	}
	std::vector<int> depth(contourCount, 0);
	for (size_t i = 0; i < contourCount; ++i) {
		int cursor = parent[i];
		while (cursor >= 0) {
			if (++depth[i] > int(contourCount))
				return false;
			cursor = parent[size_t(cursor)];
		}
	}

	std::vector<int> triangles;
	long double expectedDoubleArea = 0;
	for (size_t outerIndex = 0; outerIndex < contourCount; ++outerIndex) {
		if (depth[outerIndex] % 2 != 0)
			continue;
		std::vector<IndexedPoint2> merged = contours[outerIndex];
		if (SignedDoubleArea(merged) < 0)
			std::reverse(merged.begin(), merged.end());
		expectedDoubleArea += std::abs(areas[outerIndex]);

		std::vector<std::vector<IndexedPoint2>> holes;
		for (size_t i = 0; i < contourCount; ++i) {
			if (parent[i] != int(outerIndex))
				continue;
			std::vector<IndexedPoint2> hole = contours[i];
			if (SignedDoubleArea(hole) > 0)
				std::reverse(hole.begin(), hole.end());
			expectedDoubleArea -= std::abs(areas[i]);
			holes.push_back(std::move(hole));
		}
		std::sort(holes.begin(), holes.end(), [](const auto &a, const auto &b) {
			const auto leftmostX = [](const auto &contour) {
				double x = contour[0].point.X();
				for (const auto &point : contour)
					x = std::min(x, point.point.X());
				return x;
			};
			return leftmostX(a) < leftmostX(b);
		});
		for (const auto &hole : holes)
			if (!MergeHole(merged, hole, holes, epsilon))
				return false;
		if (!TessellateWeaklySimpleContour(merged, triangles, epsilon))
			return false;
	}

	long double triangleDoubleArea = 0;
	for (size_t i = 0; i < triangles.size(); i += 3) {
		const int ia = triangles[i];
		const int ib = triangles[i + 1];
		const int ic = triangles[i + 2];
		if (ia == ib || ib == ic || ic == ia)
			return false;
		const Point2d &a = flatPoints[size_t(ia)];
		const Point2d &b = flatPoints[size_t(ib)];
		const Point2d &c = flatPoints[size_t(ic)];
		const long double triangleArea = Orient2D(a, b, c);
		if (triangleArea <= epsilon)
			return false;
		triangleDoubleArea += triangleArea;
		const Point2d centroid(
			(a.X() + b.X() + c.X()) / 3,
			(a.Y() + b.Y() + c.Y()) / 3);
		int containingContours = 0;
		for (const auto &contour : contours)
			containingContours += PointInContour(centroid, contour, epsilon) ? 1 : 0;
		if (containingContours % 2 == 0)
			return false;
	}
	const long double areaTolerance = std::max(
		std::abs(expectedDoubleArea) * 1e-12L,
		epsilon * flatPoints.size() * 4);
	if (expectedDoubleArea <= epsilon
		|| std::abs(triangleDoubleArea - expectedDoubleArea) > areaTolerance)
		return false;

	if (delaunay)   // ear clipping guarantees validity, not quality: long contours give fans of slivers
		FlipToConstrainedDelaunay(triangles, {},
			[&](int a, int b, int c) { return Orient2D(flatPoints[size_t(a)], flatPoints[size_t(b)], flatPoints[size_t(c)]); },
			[&](int a, int b) { return (flatPoints[size_t(a)] - flatPoints[size_t(b)]).SquaredNorm(); });
	output.insert(output.end(), triangles.begin(), triangles.end());
	return true;
}

/**
 * Project and triangulate finite, coplanar 3D contours with even-odd filling.
 * Output indices address the input points flattened in contour order. Triangle
 * winding follows the Newell normal of the first non-degenerate contour.
 */
template <class CONTOUR_CONTAINER>
bool TessellatePlanarContours3(
	const CONTOUR_CONTAINER &contours,
	std::vector<int> &output,
	bool delaunay = false)
{
	if (contours.empty())
		return false;
	Point3d origin;
	bool haveOrigin = false;
	double scale = 0;
	std::vector<std::vector<Point3d>> relative;
	relative.reserve(contours.size());
	for (const auto &contour : contours) {
		if (contour.size() < 3)
			return false;
		std::vector<Point3d> relativeContour;
		relativeContour.reserve(contour.size());
		for (const auto &point : contour) {
			const Point3d absolute{
				double(point[0]), double(point[1]), double(point[2])};
			if (!std::isfinite(absolute.X()) || !std::isfinite(absolute.Y())
				|| !std::isfinite(absolute.Z()))
				return false;
			if (!haveOrigin) {
				origin = absolute;
				haveOrigin = true;
			}
			relativeContour.push_back(absolute - origin);
			scale = std::max(scale, relativeContour.back().Norm());
		}
		relative.push_back(std::move(relativeContour));
	}

	Point3d normal(0, 0, 0);
	for (const auto &contour : relative) {
		Point3d candidate(0, 0, 0);
		for (size_t i = 0; i < contour.size(); ++i)
			candidate += contour[i] ^ contour[(i + 1) % contour.size()];
		if (candidate.Norm() > scale * scale * contour.size() * 64
			* std::numeric_limits<double>::epsilon()) {
			normal = candidate;
			break;
		}
	}
	if (!(normal.Norm() > 0))
		return false;
	normal.Normalize();
	const double planeTolerance = std::max(
		scale * 1e-10, scale * 128 * std::numeric_limits<double>::epsilon());
	for (const auto &contour : relative)
		for (const Point3d &point : contour)
			if (std::abs(point * normal) > planeTolerance)
				return false;

	const Point3d absoluteNormal(
		std::abs(normal.X()), std::abs(normal.Y()), std::abs(normal.Z()));
	const int droppedAxis = absoluteNormal.X() >= absoluteNormal.Y()
			&& absoluteNormal.X() >= absoluteNormal.Z()
		? 0 : (absoluteNormal.Y() >= absoluteNormal.Z() ? 1 : 2);
	std::vector<std::vector<Point2d>> projected(relative.size());
	for (size_t i = 0; i < relative.size(); ++i) {
		projected[i].reserve(relative[i].size());
		for (const Point3d &point : relative[i]) {
			if (droppedAxis == 0)
				projected[i].push_back(Point2d(point.Y(), point.Z()));
			else if (droppedAxis == 1)
				projected[i].push_back(Point2d(point.X(), point.Z()));
			else
				projected[i].push_back(Point2d(point.X(), point.Y()));
		}
	}
	std::vector<int> triangles;
	if (!TessellatePlanarContours2(projected, triangles))
		return false;
	if (delaunay) {
		// The projection along an axis keeps orientations, so convexity is decided on it
		// exactly, but not angles: those are measured on the contours themselves.
		std::vector<Point2d> flatProjected;
		std::vector<Point3d> flatRelative;
		for (size_t i = 0; i < relative.size(); ++i) {
			flatProjected.insert(flatProjected.end(), projected[i].begin(), projected[i].end());
			flatRelative.insert(flatRelative.end(), relative[i].begin(), relative[i].end());
		}
		FlipToConstrainedDelaunay(triangles, {},
			[&](int a, int b, int c) { return planar_polygon_detail::Orient2D(flatProjected[size_t(a)], flatProjected[size_t(b)], flatProjected[size_t(c)]); },
			[&](int a, int b) { return (flatRelative[size_t(a)] - flatRelative[size_t(b)]).SquaredNorm(); });
	}
	const double projectedNormalComponent = droppedAxis == 0
		? normal.X() : (droppedAxis == 1 ? -normal.Y() : normal.Z());
	if (projectedNormalComponent < 0)
		for (size_t i = 0; i < triangles.size(); i += 3)
			std::swap(triangles[i + 1], triangles[i + 2]);
	output.insert(output.end(), triangles.begin(), triangles.end());
	return true;
}

/**
 * Triangulate one finite polygon embedded in 3D.
 *
 * The preferred path projects along the dominant component of its Newell
 * normal and uses the validated 2D ear clipper above. Real-world polygon meshes
 * frequently contain non-planar, degenerate, or self-intersecting faces, so a
 * failed projection is not a fatal error: quads are split along their shorter
 * diagonal and larger polygons use a deterministic fan. Such a fallback may
 * contain degenerate, overlapping, or crossing triangles, but always preserves
 * the input boundary order and emits exactly n-2 triangles.
 *
 * False is reserved for inputs for which triangle indices cannot be produced
 * (fewer than three points or non-finite coordinates), and leaves output
 * unchanged. Callers interested in geometry quality may inspect usedFallback.
 */
template <class POINT_CONTAINER>
bool TessellatePlanarPolygon3(
	const POINT_CONTAINER &points,
	std::vector<int> &output,
	bool *usedFallback = nullptr)
{
	if (usedFallback)
		*usedFallback = false;
	const size_t count = points.size();
	if (count < 3)
		return false;

	const Point3d origin{
		static_cast<double>(points[0][0]),
		static_cast<double>(points[0][1]),
		static_cast<double>(points[0][2])};
	std::vector<Point3d> relative(count);
	double scale = 0;
	for (size_t i = 0; i < count; ++i) {
		relative[i] = Point3d(
			double(points[i][0]) - origin.X(),
			double(points[i][1]) - origin.Y(),
			double(points[i][2]) - origin.Z());
		if (!std::isfinite(relative[i].X()) || !std::isfinite(relative[i].Y())
			|| !std::isfinite(relative[i].Z()))
			return false;
		scale = std::max(scale, relative[i].Norm());
	}
	Point3d normal(0, 0, 0);
	for (size_t i = 0; i < count; ++i)
		normal += relative[i] ^ relative[(i + 1) % count];
	const double normalNorm = normal.Norm();
	const double normalTolerance = scale * scale * count * 64
		* std::numeric_limits<double>::epsilon();
	if (normalNorm > normalTolerance) {
		normal /= normalNorm;
		std::vector<Point2d> projected;
		projected.reserve(count);
		const Point3d absoluteNormal(
			std::abs(normal.X()), std::abs(normal.Y()), std::abs(normal.Z()));
		for (size_t i = 0; i < count; ++i) {
			if (absoluteNormal.X() >= absoluteNormal.Y() && absoluteNormal.X() >= absoluteNormal.Z())
				projected.push_back(Point2d(relative[i].Y(), relative[i].Z()));
			else if (absoluteNormal.Y() >= absoluteNormal.Z())
				projected.push_back(Point2d(relative[i].X(), relative[i].Z()));
			else
				projected.push_back(Point2d(relative[i].X(), relative[i].Y()));
		}
		if (TessellatePlanarPolygon2(projected, output))
			return true;
	}

	if (usedFallback)
		*usedFallback = true;
	std::vector<int> triangles;
	triangles.reserve(3 * (count - 2));
	if (count == 4
		&& (relative[0] - relative[2]).SquaredNorm()
			> (relative[1] - relative[3]).SquaredNorm()) {
		triangles = {0, 1, 3, 1, 2, 3};
	}
	else {
		for (size_t i = 1; i + 1 < count; ++i) {
			triangles.push_back(0);
			triangles.push_back(int(i));
			triangles.push_back(int(i + 1));
		}
	}
	output.insert(output.end(), triangles.begin(), triangles.end());
	return true;
}

/*@}*/
} // end namespace
#endif
