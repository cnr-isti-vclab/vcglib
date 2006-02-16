/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2006                                                \/)\/    *
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
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.2  2006/02/15 15:40:06  corsini
Decouple SVG properties and exporter for simmetry with the other exporter

Revision 1.1  2006/02/13 16:18:09  corsini
first working version


****************************************************************************/
#ifndef __VCG_LIB_EXPORTER_SVG
#define __VCG_LIB_EXPORTER_SVG

namespace vcg 
{
	namespace edge 
	{
		namespace io 
		{


/**
 * SVG Properties.
 *
 * Support class to set the properties of the SVG exporter.
 */
class SVGProperties
{
	// definitions
public:

	//! Stroke colors.
	enum StrokeColor 
	{
		BLACK, 
		SILVER,
		GRAY,
		WHITE,
		MAROON,
		RED,
		PURPLE,
		FUCHSIA,
		GREEN,
		LIME,
		OLIVE,
		YELLOW,
		NAVY,
		BLUE,
		TEAL,
		AQUA
	};

	//! Stroke linecap types.
	enum StrokeLineCap
	{
		BUTT,
		ROUND,
		SQUARE
	};

	static const int DEFAULT_LINE_WIDTH;
	static const char * DEFAULT_LINE_COLOR;
	static const char * DEFAULT_LINE_CAP;

// private data members
private:

	//! Line width.
	int lwidth;

	//! Stroke color (see StrokeColor).
	std::string stroke_color;

	//! Stroke linecap (see StrokeLineCap).
	std::string stroke_linecap;

	//! Plane where to project the edge mesh.
	Plane3d proj;

// construction
public:

	SVGProperties()
	{
		lwidth = DEFAULT_LINE_WIDTH;
		stroke_color = DEFAULT_LINE_COLOR;
		stroke_linecap = DEFAULT_LINE_CAP;

		// default projection plane (XZ plane)
		Point3d n(0.0, 1.0, 0.0);
		proj.SetDirection(n);
		proj.SetOffset(0.0);
	}

// public methods
public:

	//! Set the line width.
	void setLineWidth(int width)
	{
		lwidth = width;
	}

	//! Set the stroke color.
	void setColor(enum StrokeColor color)
	{
		if (color == BLACK)
			stroke_color = "black";
		else if (color == SILVER)
			stroke_color = "silver";
		else if (color == GRAY)
			stroke_color = "gray";
		else if (color == WHITE)
			stroke_color = "white";
		else if (color == MAROON)
			stroke_color = "maroon";
		else if (color == RED)
			stroke_color = "red";
		else if (color == PURPLE)
			stroke_color = "purple";
		else if (color == FUCHSIA)
			stroke_color = "fuchsia";
		else if (color == GREEN)
			stroke_color = "green";
		else if (color == OLIVE)
			stroke_color = "olive";
		else if (color == LIME)
			stroke_color = "lime";
		else if (color == YELLOW)
			stroke_color = "yellow";
		else if (color == NAVY)
			stroke_color = "navy";
		else if (color == BLUE)
			stroke_color = "blue";
		else if (color == TEAL)
			stroke_color = "teal";
		else if (color == AQUA)
			stroke_color = "aqua";
	}

	//! Set the line cap style.
	void setLineCap(enum StrokeLineCap linecap)
	{
		if (linecap == BUTT)
			stroke_linecap = "butt";
		else if (linecap == ROUND)
			stroke_linecap = "round";
		else if (linecap == SQUARE)
			stroke_linecap = "square";
	}

	void setPlane(double distance, Point3d &direction)
	{
		proj.SetDirection(direction);
		proj.SetOffset(distance);
	}

// accessors
public:

	int lineWidth(){return lwidth;}
	const char * lineColor(){return stroke_color.c_str();}
	const char * lineCapStyle(){return stroke_linecap.c_str();}
	const Plane3d * projPlane(){return &proj;}

};

// DEFAULT SVG PROPERTIES
const int SVGProperties::DEFAULT_LINE_WIDTH = 2;
const char * SVGProperties::DEFAULT_LINE_COLOR = "black";
const char * SVGProperties::DEFAULT_LINE_CAP = "round";

/**
 * SVG exporter.
 *
 * This exporter save a mesh of EdgeMesh type in the SVG format.
 * Most of the features of the SVG format are not supported. 
 * The given EdgeMesh is saved as a set lines. The properties
 * of the SVG export can be set through the SVGProp class.
 */
template <class EdgeMeshType>
class ExporterSVG
{

public:

	//! Save with the default SVG properties.
	static bool Save(EdgeMeshType *mp, const char *filename)
	{
		SVGProperties properties;
		return Save(mp, filename, properties);
	}
	
	//! Save with the given SVG properties.
	static bool Save(EdgeMeshType *mp, const char *filename, SVGProperties & props)
	{
		FILE * o = fopen(filename,"w");
		if (o==NULL)
			return false;

		// initial xml tags
		fprintf(o, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
		fprintf(o, "<!-- Created with vcg library -->\n");
		fprintf(o, "<svg width=\"10cm\" height=\"10cm\" viewBox=\"0 0 1000 1000\" \n");
		fprintf(o, "  xmlns:dc=\"http://purl.org/dc/elements/1.1/\" \n");
		fprintf(o, "  xmlns:cc=\"http://web.resource.org/cc/\" \n");
		fprintf(o, "  xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" \n");
		fprintf(o, "  xmlns:svg=\"http://www.w3.org/2000/svg\" \n");
		fprintf(o, "  xmlns=\"http://www.w3.org/2000/svg\" \n");
		fprintf(o, "  id=\"svg2\"> \n");
		fprintf(o, "  <defs id=\"defs4\"/> \n");
		fprintf(o, "  <metadata id=\"metadata7\"> \n");
		fprintf(o, "    <rdf:RDF> \n");
		fprintf(o, "    <cc:Work rdf:about=\"\"> \n");
		fprintf(o, "    <dc:format>image/svg+xml</dc:format> \n");
		fprintf(o, "    <dc:type rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" /> \n");
		fprintf(o, "    </cc:Work> \n");
		fprintf(o, "    </rdf:RDF> \n");
		fprintf(o, "  </metadata> \n");

		Save(mp, o, props);

		fprintf(o, "</svg>");

		// final xml tags

		fclose(o);
		return true;
	}

	static void Save(EdgeMeshType *mp, FILE* o, SVGProperties & props)
	{
		// build vector basis (n, v1, v2)
		Point3d p1(0.0,0.0,0.0);
		Point3d p2(1.0,0.0,0.0);

		Point3d d = props.projPlane()->Direction() - p2;

		Point3d v1;
		if (d.Norm() < 0.00001)
			v1 = Point3d(0.0,0.0,1.0) - p1;
		else
			v1 = p2 - p1;

		v1.Normalize();
		Point3d v2 = v1 ^ props.projPlane()->Direction();

		std::vector<Point2f> pts;
		Point2f pmin(100000000.0f,  100000000.0f);
		Point2f pmax(-100000000.0f, -100000000.0f);

		EdgeMeshType::EdgeIterator i;
		for (i = mp->edges.begin(); i != mp->edges.end(); ++i)
		{
			Point3<EdgeMeshType::ScalarType> p1 = (*i).V(0)->P();
			Point3<EdgeMeshType::ScalarType> p2 = (*i).V(1)->P();

			Point3d p1d(p1[0], p1[1], p1[2]);
			Point3d p2d(p2[0], p2[1], p2[2]);

			// Project the line on the reference plane
			Point3d p1proj = props.projPlane()->Projection(p1d);
			Point3d p2proj = props.projPlane()->Projection(p2d);

			// Convert the 3D coordinates of the line to the uv coordinates of the plane
			Point2f pt1(static_cast<float>(p1proj * v1), static_cast<float>(p1proj * v2));
			Point2f pt2(static_cast<float>(p2proj * v1), static_cast<float>(p2proj * v2));

			pts.push_back(pt1);
			pts.push_back(pt2);

			if (pt1[0] <= pmin[0])
				pmin[0] = pt1[0];

			if (pt2[0] <= pmin[0])
				pmin[0] = pt2[0];

			if (pt1[1] <= pmin[1])
				pmin[1] = pt1[1];

			if (pt2[1] <= pmin[1])
				pmin[1] = pt2[1];

			if (pt1[0] >= pmax[0])
				pmax[0] = pt1[0];

			if (pt2[0] >= pmax[0])
				pmax[0] = pt2[0];

			if (pt1[1] >= pmax[1])
				pmax[1] = pt1[1];

			if (pt2[1] >= pmax[1])
				pmax[1] = pt2[1];
		}

		float scale = 1000.0f / std::max(pmax[0] - pmin[0], pmax[1] - pmin[1]);

		// line settings
		fprintf(o, "      <g stroke=\"%s\" stroke-linecap=\"%s\" > \n", 
			props.lineColor(), props.lineCapStyle());

		std::vector<Point2f>::iterator itPoints;
		for(itPoints = pts.begin(); itPoints != pts.end(); ++itPoints)
		{
			Point2f p1 = *itPoints;
			++itPoints;
			Point2f p2 = *itPoints;

			fprintf(o, "        <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" \n", 
				(p1[0] - pmin[0]) * scale, (p1[1] - pmin[1]) * scale,
				(p2[0] - pmin[0]) * scale, (p2[1] - pmin[1]) * scale );

			fprintf(o, "          stroke-width = \"%d\" ",props.lineWidth());
			fprintf(o, "/>\n");
		}

		fprintf(o, "  </g>\n");
	}

};


		};  // namespace io
	};  // namespace edge
};  // namespace vcg
#endif  // __VCG_LIB_EXPORTER_SVG
