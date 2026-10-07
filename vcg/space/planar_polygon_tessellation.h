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
	const auto maxCos = [&len2](int a, int b, int c) {   // larger as the smallest angle shrinks
		const double ab = len2(a, b), bc = len2(b, c), ca = len2(c, a);
		return std::max({ (ab + ca - bc) / (2 * std::sqrt(ab * ca)),
		                  (ab + bc - ca) / (2 * std::sqrt(ab * bc)),
		                  (bc + ca - ab) / (2 * std::sqrt(bc * ca)) });
	};
	const auto key = [](int a, int b) { return (std::uint64_t(std::uint32_t(a)) << 32) | std::uint32_t(b); };
	std::unordered_map<std::uint64_t, size_t> corner;   // directed edge (a, b) -> position of a in tris
	for (size_t i = 0; i < tris.size(); ++i)
		corner[key(tris[i], tris[i - i % 3 + (i + 1) % 3])] = i;
	std::deque<std::pair<int, int>> queue;
	for (size_t i = 0; i < tris.size(); ++i) {
		const int a = tris[i], b = tris[i - i % 3 + (i + 1) % 3];
		if (a < b) queue.push_back({a, b});
	}
	while (!queue.empty()) {
		const int a = queue.front().first, b = queue.front().second;
		queue.pop_front();
		if (fixed.count({std::min(a, b), std::max(a, b)}))
			continue;
		const auto e1 = corner.find(key(a, b)), e2 = corner.find(key(b, a));
		if (e1 == corner.end() || e2 == corner.end())
			continue;
		const size_t t1 = e1->second - e1->second % 3, t2 = e2->second - e2->second % 3;
		const int c = tris[t1 + (e1->second - t1 + 2) % 3], d = tris[t2 + (e2->second - t2 + 2) % 3];
		if (!(orient(c, a, d) > 0 && orient(d, b, c) > 0))
			continue;
		const double before = std::max(maxCos(a, b, c), maxCos(b, a, d));
		const double after = std::max(maxCos(c, a, d), maxCos(d, b, c));
		if (!(after < before - 1e-12))
			continue;
		for (size_t t : { t1, t2 })
			for (int k = 0; k < 3; ++k)
				corner.erase(key(tris[t + k], tris[t + (k + 1) % 3]));
		const int nt[2][3] = { { c, a, d }, { d, b, c } };
		for (int k = 0; k < 3; ++k) { tris[t1 + k] = nt[0][k]; tris[t2 + k] = nt[1][k]; }
		for (size_t t : { t1, t2 })
			for (int k = 0; k < 3; ++k)
				corner[key(tris[t + k], tris[t + (k + 1) % 3])] = t + k;
		for (auto e : { std::make_pair(a, d), std::make_pair(d, b), std::make_pair(b, c), std::make_pair(c, a) })
			queue.push_back({std::min(e.first, e.second), std::max(e.first, e.second)});
	}
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
