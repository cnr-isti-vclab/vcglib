#ifndef POLYTOQIMAGE_H
#define POLYTOQIMAGE_H
#include <QImage>
#include <QSvgGenerator>
#include <QPainter>
#include <vcg/space/point2.h>
#include <vcg/space/color4.h>
#include <vcg/space/box2.h>
#include <vcg/math/similarity2.h>

///this class is used to pass global 
///parameters to the polygonal dumper
class PolyDumperParam
{
public:
	///the backgrround color
	vcg::Color4b backgroundColor;
	///true if the polygons must be filled with color
	bool fill;
	///true if the filling color is random
	bool randomColor;
	///the filling color of polygons, used only if randomColor==false
	vcg::Color4b FillColor;
	
	///dimension of the image (in PNG are pixels, while in SV is the workspace in points)
	int width;
	int height;
	///DPi resolution, used only for SVG
	int dpi;
	
	///default contructor
	PolyDumperParam()
	{
		backgroundColor = vcg::Color4b::Gray;
		width=1024;
		height=1024;
		dpi=72;
		fill=false;
		randomColor=true;
		FillColor=vcg::Color4b(0,0,0,255);
	}
};

///this class is used to draw polygons on an image could be vectorial or not
class PolyDumper
{
	///this class draw a black mask fora given polygon, cenetered and scaled to fit with 
	///the image size, it return the transformation to tranform back the polygon to 2D space
	static void DrawPolygonMask(const std::vector< std::vector<vcg::Point2f> > &polyVec,QImage &img,
								vcg::Similarity2f &ret,const vcg::Similarity2f &trans);
	
	///return the max radius of a point inside a polygon ,given the mask image
	///actually it evaluate the maximum bounding box
	static int getMaxMaskRadius(int x,int y,QImage &img);
	
	///return the point inside the polygon with the bigger distance to the border,
	///this is used to write labels within the polygon, it handle polygons with holes too
	static vcg::Point2f GetIncenter(const std::vector< std::vector<vcg::Point2f> > &polyVec,
										const vcg::Similarity2f &tra1,int &radius,int resolution=100);

	static void rectSetToPolySet(std::vector< vcg::Box2f > &rectVec, std::vector< std::vector<vcg::Point2f> > &polyVec);

public:
	///write a polygon on a PNG file, format of the polygon is vector of vector of contours...nested contours are holes
	///takes the name of the image in input, the set of polygons, the set of per polygons transformation, 
	///the label to be written and the global parameter for drawing style
	static void dumpPolySetPNG(const char * imageName, std::vector< std::vector< std::vector<vcg::Point2f> > > &polyVecVec,
							std::vector<vcg::Similarity2f> &trVec, std::vector<std::string> &labelVec, PolyDumperParam &pp);
	//write a polygon on a SVG file, format of the polygon is vector of vector of contours...nested contours are holes
	///takes the name of the image in input, the set of polygons, the set of per polygons transformation, 
	///the label to be written and the global parameter for drawing style
	static void dumpPolySetSVG(const char * imageName, std::vector< std::vector< std::vector<vcg::Point2f> > > &polyVecVec,
							std::vector<vcg::Similarity2f> &trVec, std::vector<std::string> &labelVec, PolyDumperParam &pp);
	static void dumpPolySetPNG(const char * imageName, std::vector< std::vector< std::vector<vcg::Point2f> > > &polyVecVec,
							std::vector<vcg::Similarity2f> &trVec, PolyDumperParam &pp);
	static void dumpPolySetSVG(const char * imageName, std::vector< std::vector< std::vector<vcg::Point2f> > > &polyVecVec, 
							std::vector<vcg::Similarity2f> &trVec, PolyDumperParam &pp);
	static void dumpPolySetPNG(const char * imageName, std::vector<  std::vector<vcg::Point2f> > &polyVecVec, 
							std::vector<vcg::Similarity2f> &trVec, PolyDumperParam &pp);
	static void dumpPolySetSVG(const char * imageName, std::vector<  std::vector<vcg::Point2f> > &polyVecVec, 
							std::vector<vcg::Similarity2f> &trVec, PolyDumperParam &pp);
};
#endif // POLYTOQIMAGE_H
