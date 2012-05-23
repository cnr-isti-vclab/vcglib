

#include "PolyToQImage.h"
#include <wrap/qt/col_qt_convert.h>
using namespace vcg;
using namespace std;

void PolyDumper::rectSetToPolySet(vector< Box2f > &rectVec, vector< vector<Point2f> > &polyVec)
{
	polyVec.clear();
	for(size_t i=0;i<rectVec.size();++i)
	{
		Box2f &b=rectVec[i];
		polyVec.resize(polyVec.size()+1);
		polyVec.back().push_back(b.min);
		polyVec.back().push_back(Point2f(b.max[0],b.min[1]));
		polyVec.back().push_back(b.max);
		polyVec.back().push_back(Point2f(b.min[0],b.max[1]));
	}
}


void PolyDumper::dumpPolySetSVG(const char * imageName, vector< vector<Point2f> > &polyVec, vector<Similarity2f> &trVec, PolyDumperParam &pp)
{
	vector< vector< vector<Point2f> > > polyVecVec(polyVec.size());
	for(size_t i=0;i<polyVec.size();++i)
	{
		polyVecVec[i].resize(1);
		polyVecVec[i][0]=polyVec[i];
	}
	dumpPolySetSVG(imageName,polyVecVec,trVec,pp);
}

void PolyDumper::dumpPolySetPNG(const char * imageName, vector< vector<Point2f> > &polyVec, vector<Similarity2f> &trVec, PolyDumperParam &pp)
{
	vector< vector< vector<Point2f> > > polyVecVec(polyVec.size());
	for(size_t i=0;i<polyVec.size();++i)
	{
		polyVecVec[i].resize(1);
		polyVecVec[i][0]=polyVec[i];
	}
	dumpPolySetPNG(imageName,polyVecVec,trVec,pp);
}

void PolyDumper::dumpPolySetPNG(const char * imageName, vector< vector< vector<Point2f> > > &polyVecVec, vector<Similarity2f> &trVec, PolyDumperParam &pp)
{
	vector<string> labelVec(polyVecVec.size());
	dumpPolySetPNG(imageName,polyVecVec,trVec,labelVec,pp);
}

void PolyDumper::dumpPolySetSVG(const char * imageName, vector< vector< vector<Point2f> > > &polyVecVec, vector<Similarity2f> &trVec, PolyDumperParam &pp)
{
	vector<string> labelVec(polyVecVec.size());
	dumpPolySetSVG(imageName,polyVecVec,trVec,labelVec,pp);
}

///this class draw a black mask fora given polygon, cenetered and scaled to fit with 
///the image size, it return the transformation to tranform back the polygon to 2D space
void PolyDumper::DrawPolygonMask(const vector< vector<Point2f> > &polyVec,
								QImage &img,
								Similarity2f &ret,
								const Similarity2f &trans)
{

	///RASTERIZE THE POLYGON
	QPainter painter;
	painter.begin(&img);
	QBrush bb;
	bb.setStyle(Qt::SolidPattern);

	int resolution=img.width();

	QPen qp;
	qp.setWidthF(0);
	///find the BB
	vcg::Box2f bbox;
	for(size_t i=0;i<polyVec.size();++i)
		for(size_t j=0;j<polyVec[i].size();++j)
		{
			Point2f pp=polyVec[i][j];
			pp.Rotate(trans.rotRad);
			bbox.Add(pp);
		}
		vcg::Point2f pos=bbox.Center();
		float scaleFact=(resolution/bbox.Diag());

		///SET THE TRANSFORMATION TO CENTER WRT BBOX
		painter.resetTransform();
		painter.translate(resolution/2,resolution/2);
		painter.rotate(math::ToDeg(trans.rotRad));
		painter.scale(scaleFact,scaleFact);
		painter.translate(-pos.V(0),-pos.V(1));

		///SET TO RETURN BACK
		ret.rotRad=0;
		ret.sca=scaleFact;
		ret.tra[0]=-pos.V(0)*scaleFact + resolution/2;
		ret.tra[1]=-pos.V(1)*scaleFact + resolution/2;

		///DRAW THE POLYGON
		QPainterPath QPP;
		for(size_t i=0;i<polyVec.size();++i)
		{
			QVector<QPointF> ppQ;
			for(size_t j=0;j<polyVec[i].size();++j)
			{
				Point2f pp=polyVec[i][j];
				//pp.Rotate(trans.rotRad);
				ppQ.push_back(QPointF(pp[0],pp[1]));
			}
			ppQ.push_back(QPointF(polyVec[i][0][0],polyVec[i][0][1]));
			QPP.addPolygon(QPolygonF(ppQ));
		}

		painter.setBrush(bb);
		painter.setPen(qp);
		painter.drawPath(QPP);
		painter.end();
		//img.save("test.png");
}

///return the max radius of a point inside a polygon ,given the mask image
///actually it evaluate the maximum bounding box
int PolyDumper::getMaxMaskRadius(int x,int y,QImage &img)
{
	int sizeY=img.size().height();
	int sizeX=img.size().width();
	int Max_radiusX=std::min(abs(x-sizeX),x);
	int Max_radiusY=std::min(abs(y-sizeY),y);
	int Max_radius=std::min(Max_radiusX,Max_radiusY);
	int val=qGray(img.pixel(x,y));
	///if is outside
	assert(val!=255);

	int curr_radius=1;
	while (curr_radius<Max_radius)
	{
		vcg::Box2i currBB=vcg::Box2i(vcg::Point2i(x,y),vcg::Point2i(x,y));
		vcg::Box2i InsideBB=currBB;
		currBB.Offset(curr_radius);
		InsideBB.Offset(curr_radius-1);

		///check on borders
		for (int x0=currBB.min.X();x0<=currBB.max.X();x0++)
			for (int y0=currBB.min.Y();y0<=currBB.max.Y();y0++)
			{
				if (InsideBB.IsIn(vcg::Point2i(x0,y0)))continue;
				int val=qGray(img.pixel(x0,y0));
				///if is outside
				if (val==255)
					return curr_radius;
			}
			curr_radius++;
	}
	return curr_radius;
}

///return the point inside the polygon with the bigger distance to the border,
///this is used to write labels within the polygon, it handle polygons with holes too
vcg::Point2f PolyDumper::GetIncenter(const vector< vector<Point2f> > &polyVec,
											const Similarity2f &tra1,
											int &radius,
											int resolution)
{
	///INITIALIZE THE IMAGE
	QImage img(resolution,resolution,QImage::Format_RGB32);
	img.fill(vcg::ColorConverter::ToQColor(vcg::Color4b(255,255,255,255)).rgb());
	Similarity2f tra0;
	///DRAW THE MASK
	DrawPolygonMask(polyVec,img,tra0,tra1);
	//img = img.mirrored(false,true);

	///THEN FIND THE CENTROID
	float Maxradius=-1;
	int sizeY=img.size().height();
	int sizeX=img.size().width();
	vcg::Point2i incenter=vcg::Point2i(sizeX/2,sizeY/2);
	///THEN GO OVER ALL VERTICES
	for (int x=0;x<sizeX;x++)
		for (int y=0;y<sizeY;y++)
		{
			///flipped because of 
			///same artifacts of Mquads
			int val=qGray(img.pixel(x,y));
			///if is outside
			if (val==255)continue;
			int curr_radius=getMaxMaskRadius(x,y,img);
			if (curr_radius>Maxradius)
			{
				Maxradius=curr_radius;
				incenter=vcg::Point2i(x,y);
			}
		}

		vcg::Point2f incenterf;
		incenterf.Import(incenter);
//		/FIRST TRASNFORMATION
//		QPainter painter;
//		painter.begin(&img);
//		//painter.fillRect(incenter.V(0)-Maxradius,incenter.V(1)-Maxradius,Maxradius*2,Maxradius*2,Qt::red);
////		painter.fillRect(incenter.V(0)-1,incenter.V(1)-1,2,2,Qt::red);
//		//painter.drawPoint(incenter.V(0),incenter.V(1));
//		painter.end();
//		static int num=0;
//		num++;
//		char path[100];
//		sprintf(path,"mask%d.png",num);
//		img.save(path);


		incenterf.X()-=tra0.tra[0];
		incenterf.Y()-=tra0.tra[1];
		incenterf*=1.0/tra0.sca;
		incenterf*=tra1.sca;

		///SECOND TRANSFORMATION
		assert(Maxradius>0);
		radius=Maxradius;
		return incenterf;
}


///write a polygon on a PNG file, format of the polygon is vector of vector of contours...nested contours are holes
///takes the name of the image in input, the set of polygons, the set of per polygons transformation, 
///the label to be written and the global parameter for drawing style
void  PolyDumper::dumpPolySetPNG(const char * imageName, 
										vector< vector< vector<Point2f> > > &polyVecVec, 
										vector<Similarity2f> &trVec, 
										vector<string> &labelVec, 
										PolyDumperParam &pp)
{
	///SET THE FONT
	assert(polyVecVec.size() == trVec.size());
	int fontSize=ceil(std::max(pp.width,pp.height)/100.0);
	QFont qf("courier",fontSize);

	///SET THE DRAWING SIZE
	QImage img(pp.width,pp.height,QImage::Format_RGB32);
	img.fill(vcg::ColorConverter::ToQColor( pp.backgroundColor).rgb());

	///SETUP OF DRAWING PROCEDURE
	QPainter painter;
	painter.begin(&img);
	QBrush bb;
	if (pp.fill)
		bb.setStyle(Qt::SolidPattern);
	QPen qp;
	qp.setWidthF(0);

	for(size_t i=0;i<polyVecVec.size();++i)
	{
		///SET THE CURRENT TRANSFORMATION
		painter.resetTransform();
		painter.translate(trVec[i].tra[0],trVec[i].tra[1]);
		painter.rotate(math::ToDeg(trVec[i].rotRad));
		painter.scale(trVec[i].sca,trVec[i].sca);
		QPainterPath QPP;

		for(size_t jj=0;jj<polyVecVec[i].size();++jj)
		{
			QVector<QPointF> ppQ;
			for(size_t j=0;j<polyVecVec[i][jj].size();++j)
			{
				Point2f pp=polyVecVec[i][jj][j];
				ppQ.push_back(QPointF(pp[0],pp[1]));
			}
			ppQ.push_back(QPointF(polyVecVec[i][jj][0][0],polyVecVec[i][jj][0][1]));
			QPP.addPolygon(QPolygonF(ppQ));
		}
		///FIND THE BARYCENTER
		int radius;
		Point2f bc;
		bc=GetIncenter(polyVecVec[i],trVec[i],radius,100);

		if (pp.randomColor)
			bb.setColor(vcg::ColorConverter::ToQColor(Color4b::Scatter(polyVecVec.size(),i)));
		else
			bb.setColor(vcg::ColorConverter::ToQColor(pp.FillColor));

		painter.setBrush(bb);
		painter.setPen(qp);
		painter.drawPath(QPP);

		///DRAW THE TEXT
		painter.setFont(qf);
		painter.resetTransform();
		painter.translate(trVec[i].tra[0],trVec[i].tra[1]);

		painter.drawText(bc[0]-radius,bc[1]-radius,radius*2,radius*2,Qt::AlignHCenter|Qt::AlignCenter,QString(labelVec[i].c_str()));
	}
	painter.end();
	img.save(imageName);
}


///write a polygon on a SVG file, format of the polygon is vector of vector of contours...nested contours are holes
///takes the name of the image in input, the set of polygons, the set of per polygons transformation, 
///the label to be written and the global parameter for drawing style
void PolyDumper::dumpPolySetSVG(const char * imageName, 
									   vector< vector< vector<Point2f> > > &polyVecVec, 
									   vector<Similarity2f> &trVec, 
									   vector< vector< string> > &labelVecVec,
									   vector<vector<Point2f> > &labelPosVecVec,
									   vector<vector<float> >&labelRadVecVec,
									   PolyDumperParam &pp)
{
	assert(polyVecVec.size() == trVec.size());


	///SET THE FONT
	int fontSize;
	if(pp.fontSize==0) fontSize=ceil(std::max(pp.width,pp.height)/200.0);
	else fontSize=pp.fontSize;
	QFont qf("courier",fontSize);
	QSvgGenerator svg;
	svg.setFileName(imageName);

	///SET THE DRAWING SIZE
	svg.setSize(QSize(pp.width,pp.height));
	svg.setViewBox(QRect(0, 0, pp.width,pp.height));
	svg.setResolution(int(pp.dpi));//

	///SETUP OF DRAWING PROCEDURE
	QPainter painter;
	painter.begin(&svg);
	QBrush bb;
	if (pp.fill)
		bb.setStyle(Qt::SolidPattern);

	QPen qp;

	///SET THE GLOBAL SCALING FACTOR
	float Scalesvg=1.f/(float)trVec[0].sca;
	qp.setWidthF(Scalesvg);
	for(size_t i=0;i<polyVecVec.size();++i)
	{
		///SET THE CURRENT TRANSFORMATION
		painter.resetTransform();
		painter.translate(trVec[i].tra[0],trVec[i].tra[1]);
		painter.rotate(math::ToDeg(trVec[i].rotRad));
		painter.scale(trVec[i].sca,trVec[i].sca);
		QPainterPath QPP;

		for(size_t jj=0;jj<polyVecVec[i].size();++jj)
		{
			QVector<QPointF> ppQ;
			for(size_t j=0;j<polyVecVec[i][jj].size();++j)
			{
				Point2f pp=polyVecVec[i][jj][j];
				ppQ.push_back(QPointF(pp[0],pp[1]));
			}
			ppQ.push_back(QPointF(polyVecVec[i][jj][0][0],polyVecVec[i][jj][0][1]));
			QPP.addPolygon(QPolygonF(ppQ));
		}
		///FIND THE INCENTER

		if (pp.randomColor)
			bb.setColor(vcg::ColorConverter::ToQColor(Color4b::Scatter(polyVecVec.size(),i)));
		else
			bb.setColor(vcg::ColorConverter::ToQColor(pp.FillColor));

		///DRAW THE POLYGON
		painter.setBrush(bb);
		painter.setPen(qp);
		painter.drawPath(QPP);

		///DRAW THE TEXT
		painter.setFont(qf);
		float radius;
		int radiusInt;
		Point2f bc;
		// if we do not have labelPos use the old method of empty disk.
		if(labelPosVecVec.empty()) bc=GetIncenter(polyVecVec[i],trVec[i],radiusInt);
		radius = radiusInt;
		for(size_t labelInd=0;labelInd<labelVecVec[i].size();++labelInd)
		{
		  if(!labelPosVecVec.empty())
		  {
			bc= labelPosVecVec[i][labelInd];
			bc.Rotate(trVec[i].rotRad);
			bc *= trVec[i].sca;
			radius=labelRadVecVec[i][labelInd];
			radius *= trVec[i].sca;
		  }
		  painter.resetTransform();
		  painter.translate(trVec[i].tra[0],trVec[i].tra[1]);
		  painter.drawText(bc[0]-radius,bc[1]-radius,radius*2,radius*2,Qt::AlignHCenter|Qt::AlignCenter,QString(labelVecVec[i][labelInd].c_str()));
		}
	}
	painter.end();
}

void PolyDumper::dumpPolySetSVG(const char * imageName,
									   vector< vector< vector<Point2f> > > &polyVecVec,
									   vector<Similarity2f> &trVec,
									   vector< string > &labelVec,
									   PolyDumperParam &pp)
{
  vector< vector< string> > labelVecVec(labelVec.size());
  vector< vector<Point2f> > labelPosVec;
  vector< vector<float> >labelRadVec;
  for(size_t i=0;i<labelVec.size();++i)
  {
    labelVecVec[i].push_back(labelVec[i]);
  }
  dumpPolySetSVG(imageName,polyVecVec,trVec,labelVecVec,labelPosVec,labelRadVec,pp);
}

