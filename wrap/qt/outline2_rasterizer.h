#ifndef QTPOLYRASTERIZER_H
#define QTPOLYRASTERIZER_H

#include <QImage>
#include <QSvgGenerator>
#include <QPainter>
#include <vcg/space/point2.h>
#include <vcg/space/color4.h>
#include <vcg/space/box2.h>
#include <vcg/math/similarity2.h>

///this class is used to pass global
///parameters to the polygonal dumper
namespace vcg{
class RasterizedOutline2;
}

///this class is used to draw polygons on an image could be vectorial or not
class QtPolyRasterizer
{
public:
    static void rasterize(vcg::RasterizedOutline2 &poly,
                          float scaleFactor,
                          int rast_i, int rotationNum, int cellSize);

    static std::vector<std::vector<int> > rotateGridCWise(std::vector< std::vector<int> >& inGrid);

};
#endif // QTPOLYRASTERIZER_H
