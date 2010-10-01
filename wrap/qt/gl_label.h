/****************************************************************************
* MeshLab                                                           o o     *
* An extendible mesh processor                                    o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005, 2009                                          \/)\/    *
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
#ifndef VCG_WRAP_QT_GLLABEL_H
#define VCG_WRAP_QT_GLLABEL_H

#include <QMessageBox>
#include <wrap/qt/col_qt_convert.h>
#include <wrap/qt/checkGLError.h>
namespace vcg
{

  class glLabel
  {
	public:


    static void render(QPainter *painter, const vcg::Point3f &p, const QString &text)
    {
      QFont qFont;
      qFont.setStyleStrategy(QFont::NoAntialias);
      qFont.setFamily("Helvetica");
      qFont.setPixelSize(12);
      vcg::Color4b color(vcg::Color4b::White);
      render(painter,p,text,qFont,color);
    }

    static void render(QPainter *painter, const vcg::Point3f &p, const QString &text, QFont qFont, vcg::Color4b color, float angle =0,bool rightAlign=false)
    {
      GLdouble model[16];
      GLdouble proj[16];
      GLint view[4];

      glGetDoublev(GL_MODELVIEW_MATRIX, model);
      glGetDoublev(GL_PROJECTION_MATRIX, proj);
      glGetIntegerv(GL_VIEWPORT, view);
      GLdouble winx,winy,winz;

      gluProject(p[0],p[1],p[2],model,proj,view,&winx,&winy,&winz);

      glPushAttrib(GL_ENABLE_BIT);
      glDisable(GL_DEPTH_TEST);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      QFontMetrics qfm(qFont);
      QRect textBox = qfm.boundingRect(text);
      painter->endNativePainting();
      painter->save();
      painter->setRenderHint(QPainter::TextAntialiasing);
      painter->setPen(vcg::ColorConverter::ToQColor(color));
      painter->setFont(qFont);
      //qDebug("Printing at %f %f (%f %f %f) '%s'",winx,winy,p[0],p[1],p[2],qPrintable(text));
      painter->translate(QPointF(winx,view[3]-winy));
      painter->rotate(angle);
      //painter->drawText(QPointF(0,0),text);
      //painter->drawText(QRect(),Qt::AlignLeft+Qt::AlignVCenter,text);
      //painter->drawText(QRect(view[0],view[1],view[2],view[3]),Qt::AlignLeft+Qt::AlignVCenter,text);
      //QPoint base(0,textBox.height()/2);
      QPoint base(0,qfm.ascent()/2);
      if(rightAlign)
        base.setX(-textBox.width() -qfm.maxWidth());
      painter->drawText(base,text);
      //painter->drawText(QPointF(winx,view[3]-winy),text);
      painter->restore();
      painter->beginNativePainting();
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
      glPopAttrib();
      checkGLError::qDebug("myRenderText");
    }



    static void render(QPainter *painter, const vcg::Point3d &p, const QString &text)
    { render(painter,Point3f::Construct(p),text); }
    static void render(QPainter *painter, const vcg::Point3d &p, const QString &text, QFont qFont, vcg::Color4b color)
    { render(painter,Point3f::Construct(p),text,qFont,color); }

  };
} // end namespace

#endif
