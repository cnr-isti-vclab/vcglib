#include <wrap/qt/outline2_rasterizer.h>
#include <wrap/qt/col_qt_convert.h>
#include <vcg/space/color4.h>
#include <wrap/qt/col_qt_convert.h>

#include <fstream>

using namespace vcg;
using namespace std;

void QtOutline2Rasterizer::rasterize(RasterizedOutline2 &poly,
                                 float scale,
                                 int rast_i,
                                 int rotationNum,
                                 int gutterWidth)
{

    gutterWidth *= 2; // since the brush is centered on the outline multiply the given value by 2

    float rotRad = M_PI*2.0f*float(rast_i) / float(rotationNum);

    //get polygon's BB, rotated according to the input parameter
    Box2f bb;
    vector<Point2f> pointvec = poly.getPoints();
    for(size_t i=0;i<pointvec.size();++i) {
        Point2f pp=pointvec[i];
        pp.Rotate(rotRad);
        bb.Add(pp);
    }

    //create the polygon to print it
    QVector<QPointF> points;
    vector<Point2f> newpoints = poly.getPoints();
    for (size_t i = 0; i < newpoints.size(); i++) {
        points.push_back(QPointF(newpoints[i].X(), newpoints[i].Y()));
    }

    // Compute the raster space size by rounding up the scaled bounding box size
    // and adding the gutter width.
    int sizeX = (int)ceil(bb.DimX()*scale);
    int sizeY = (int)ceil(bb.DimY()*scale);
    int safetyBuffer = 2;
    sizeX += (gutterWidth + safetyBuffer);
    sizeY += (gutterWidth + safetyBuffer);

    QImage img(sizeX,sizeY,QImage::Format_RGB32);
    QColor backgroundColor(Qt::transparent);
    img.fill(backgroundColor);

    ///SETUP OF DRAWING PROCEDURE
    QPainter painter;
    painter.begin(&img);
    {
        QBrush br;
        br.setStyle(Qt::SolidPattern);
        br.setColor(Qt::yellow);

        QPen qp;
        qp.setWidthF(0);
        qp.setWidth(gutterWidth);
        qp.setCosmetic(true);
        qp.setColor(Qt::yellow);
        qp.setJoinStyle(Qt::MiterJoin);
        qp.setMiterLimit(0);

        painter.setBrush(br);
        painter.setPen(qp);

        painter.resetTransform();
        painter.translate(QPointF(-(bb.min.X()*scale) + (gutterWidth + safetyBuffer)/2.0f, -(bb.min.Y()*scale) + (gutterWidth + safetyBuffer)/2.0f));
        painter.rotate(math::ToDeg(rotRad));
        painter.scale(scale,scale);

        painter.drawPolygon(QPolygonF(points));
    }
    painter.end();

    // workaround/hack to avoid ``disappearing'' primitives: use a cosmetic pen to
    // draw the poly boundary.
    // The proper way to do this would be to use conservative reasterization, which
    // Qt doesn't seem to support
    std::vector<QPointF> lines;
    for (int i = 1; i < points.size(); ++i) {
        lines.push_back(points[i-1]);
        lines.push_back(points[i]);
    }
    lines.push_back(points.back());
    lines.push_back(points.front());

    painter.begin(&img);
    {
        QBrush br;
        br.setStyle(Qt::SolidPattern);
        br.setColor(Qt::yellow);

        QPen qp;
        qp.setWidthF(0);
        qp.setWidth(std::max(1, gutterWidth));
        qp.setCosmetic(true);
        qp.setColor(Qt::yellow);

        painter.setBrush(br);
        painter.setPen(qp);

        painter.resetTransform();
        painter.translate(QPointF(-(bb.min.X()*scale) + (gutterWidth + safetyBuffer)/2.0f, -(bb.min.Y()*scale) + (gutterWidth + safetyBuffer)/2.0f));
        painter.rotate(math::ToDeg(rotRad));
        painter.scale(scale,scale);

        //painter.drawPoints(QPolygonF(points));
        painter.drawLines(lines.data(), lines.size()/2);
    }
    painter.end();

    // Cropping

    /*
    // Slower version
    int minX = img.width();
    int minY = img.height();
    int maxX = -1;
    int maxY = -1;

    for (int i = 0; i < img.height(); ++i) {
        const QRgb *line = reinterpret_cast<const QRgb*>(img.scanLine(i));
        for (int j = 0; j < img.width(); ++j) {
            if (line[j] != backgroundColor.rgb()) {
                if (j < minX) minX = j;
                if (j > maxX) maxX = j;
                if (i < minY) minY = i;
                if (i > maxY) maxY = i;
            }
        }
    }
    */

    int minX = img.width();
    int minY = img.height();
    int maxX = 0;
    int maxY = 0;

    for (int i = 0; i < img.height(); ++i) {
        const QRgb *line = reinterpret_cast<const QRgb*>(img.scanLine(i));
        for (int j = 0; j < img.width(); ++j) {
            if (line[j] != backgroundColor.rgb()) {
                minY = i;
                break;
            }
        }
        if (minY < img.height()) break;
    }

    for (int i = img.height() - 1; i >= 0; --i) {
        const QRgb *line = reinterpret_cast<const QRgb*>(img.scanLine(i));
        for (int j = 0; j < img.width(); ++j) {
            if (line[j] != backgroundColor.rgb()) {
                maxY = i;
                break;
            }
        }
        if (maxY > 0) break;
    }

    for (int i = minY; i <= maxY; ++i) {
        const QRgb *line = reinterpret_cast<const QRgb*>(img.scanLine(i));
        for (int j = 0; j < minX; ++j)
            if (line[j] != backgroundColor.rgb() && j < minX) {
                minX = j;
                break;
            }
        for (int j = img.width() - 1; j >= maxX; --j)
            if (line[j] != backgroundColor.rgb() && j > maxX) {
                maxX = j;
                break;
            }
    }

    assert (minX <= maxX && minY <= maxY);

    int imgW = (maxX - minX) + 1;
    int imgH = (maxY - minY) + 1;

    {
        QImage imgcp = img.copy(0, 0, img.width(), img.height());
        img = imgcp.copy(minX, minY, imgW, imgH);
    }

    //create the first grid, which will then be rotated 3 times.
    //we will reuse this grid to create the rasterizations corresponding to this one rotated by 90/180/270°
    vector<vector<int> > tetrisGrid;
    QRgb yellow = QColor(Qt::yellow).rgb();
    tetrisGrid.resize(img.height());
    for (int k = 0; k < img.height(); k++) {
        tetrisGrid[k].resize(img.width(), 0);
    }
    for (int y = 0; y < img.height(); y++) {
        const uchar* line = img.scanLine(y);
        for(int x = 0; x < img.width(); ++x) {
            if (((QRgb*)line)[x] == yellow) {
                tetrisGrid[y][x] = 1;
            }
        }
    }

    //create the 4 rasterizations (one every 90°) using the discrete representation grid we've just created
    int rotationOffset = rotationNum/4;
    for (int j = 0; j < 4; j++) {
        if (j != 0)  {
            tetrisGrid = rotateGridCWise(tetrisGrid);
        }
        //add the grid to the poly's vector of grids
        poly.getGrids(rast_i + rotationOffset*j) = tetrisGrid;

        //initializes bottom/left/deltaX/deltaY vectors of the poly, for the current rasterization
        poly.initFromGrid(rast_i + rotationOffset*j);
    }
}

// rotates the grid 90 degree clockwise (by simple swap)
// used to lower the cost of rasterization.
vector<vector<int> > QtOutline2Rasterizer::rotateGridCWise(vector< vector<int> >& inGrid) {
    vector<vector<int> > outGrid(inGrid[0].size());
    for (size_t i = 0; i < inGrid[0].size(); i++) {
        outGrid[i].reserve(inGrid.size());
        for (size_t j = 0; j < inGrid.size(); j++) {
            outGrid[i].push_back(inGrid[inGrid.size() - j - 1][i]);
        }
    }
    return outGrid;
}
