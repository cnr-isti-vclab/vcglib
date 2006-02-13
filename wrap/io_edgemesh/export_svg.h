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
 * SVG exporter.
 *
 * This exporter save a mesh of EdgeMesh type in the SVG format.
 * Most of the features of the SVG format are not supported. 
 * The given EdgeMesh is saved as a set lines. It is possible to 
 * set the width and the color of the lines.
 */
template <class EdgeMeshType>
class ExporterSVG
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

// private data members
private:

	//! Line width.
	int lwidth;

	//! Stroke color (see StrokeColor).
	std::string stroke_color;

	//! Stroke linecap (see StrokeLineCap).
	std::string stroke_linecap;

// public methods
public:

	ExporterSVG(void)
	{
		lwidth = 5;
		stroke_color = "black";
		stroke_linecap = "round";
	}

	static void setLineWidth(int width)
	{
		lwidth = width;
	}

	static void setColor(enum StrokeColor color)
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

	static void setLineCap(enum StrokeLineCap linecap)
	{
		if (linecap == BUTT)
			stroke_linecap = "butt";
		else if (linecap == ROUND)
			stroke_linecap = "round";
		else if (linecap == SQUARE)
			stroke_linecap = "square";
	}
	
	bool Save(EdgeMeshType  *mp, const char * filename)
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

		Save(mp,o);

		fprintf(o, "</svg>");

		// final xml tags

		fclose(o);
		return true;
	}

	void Save(EdgeMeshType *mp, FILE* o)
	{
		EdgeMeshType::EdgeIterator i;

		Point3f pmin = mp->bbox.min;
		float scale = 1000.0f / max(mp->bbox.DimX(), mp->bbox.DimZ());

		// line settings
		fprintf(o, "      <g stroke=\"%s\" stroke-linecap=\"%s\" > \n", 
			stroke_color.c_str(), stroke_linecap.c_str());

		for(i = mp->edges.begin(); i != mp->edges.end(); ++i)
		{
			Point3f p1 = (*i).V(0)->P();
			Point3f p2 = (*i).V(1)->P();

			fprintf(o, "        <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" \n", 
				(p1[0] - pmin[0]) * scale, (p1[2] - pmin[2]) * scale,
				(p2[0] - pmin[0]) * scale, (p2[2] - pmin[2]) * scale );

			fprintf(o, "          stroke-width = \"%d\" ",lwidth);
			fprintf(o, "/>\n");
		}

		fprintf(o, "  </g>\n");
	}

};


		};  // namespace io
	};  // namespace edge
};  // namespace vcg
#endif  // __VCG_LIB_EXPORTER_SVG
