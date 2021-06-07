#ifndef SHOT_QT_H
#define SHOT_QT_H

#include <QDomDocument>
#include <QStringList>

#include <vcg/math/shot.h>

/** This function read a shot from a parsed XML node.

  */

template <class ShotType>
    bool ReadShotFromQDomNode(
        ShotType &shot, /// the shot that will contain the read node
        const QDomNode &node) /// The XML node to be read
{
  typedef typename ShotType::ScalarType ScalarType;
  typedef vcg::Point3<typename ShotType::ScalarType> Point3x;
  if(QString::compare(node.nodeName(),"VCGCamera")==0)
  {
    QDomNamedNodeMap attr = node.attributes();
    if (attr.contains("BinaryData") && attr.namedItem("BinaryData").nodeValue().toInt() == 1)
    {
      Point3x tra;
      QString str = attr.namedItem("TranslationVector").nodeValue();
      QByteArray value = QByteArray::fromBase64(str.toLocal8Bit());
      memcpy(tra.V(), value.data(), sizeof(ScalarType) * 3);
      shot.Extrinsics.SetTra(-tra);

      vcg::Matrix44<ScalarType> rot;
      str = attr.namedItem("RotationMatrix").nodeValue();
      value = QByteArray::fromBase64(str.toLocal8Bit());
      memcpy(rot.V(), value.data(), sizeof(ScalarType) * 16);
      shot.Extrinsics.SetRot(rot);

      vcg::Camera<ScalarType> &cam = shot.Intrinsics;
      if (attr.contains("CameraType")) cam.cameraType = attr.namedItem("CameraType").nodeValue().toInt();
      memcpy(&cam.FocalMm, QByteArray::fromBase64(attr.namedItem("FocalMm").nodeValue().toLocal8Bit()).data(), sizeof(ScalarType));
      memcpy(&cam.ViewportPx, QByteArray::fromBase64(attr.namedItem("ViewportPx").nodeValue().toLocal8Bit()).data(), sizeof(int) * 2);
      memcpy(&cam.CenterPx, QByteArray::fromBase64(attr.namedItem("CenterPx").nodeValue().toLocal8Bit()).data(), sizeof(ScalarType)*2);
      memcpy(&cam.PixelSizeMm, QByteArray::fromBase64(attr.namedItem("PixelSizeMm").nodeValue().toLocal8Bit()).data(), sizeof(ScalarType) * 2);
      memcpy(&cam.k, QByteArray::fromBase64(attr.namedItem("LensDistortion").nodeValue().toLocal8Bit()).data(), sizeof(ScalarType) * 2);
    }
    else
    {
      Point3x tra;
      tra[0] = attr.namedItem("TranslationVector").nodeValue().section(' ', 0, 0).toDouble();
      tra[1] = attr.namedItem("TranslationVector").nodeValue().section(' ', 1, 1).toDouble();
      tra[2] = attr.namedItem("TranslationVector").nodeValue().section(' ', 2, 2).toDouble();
      shot.Extrinsics.SetTra(-tra);

      vcg::Matrix44<ScalarType> rot;
      QStringList values = attr.namedItem("RotationMatrix").nodeValue().split(" ", QString::SkipEmptyParts);
      for (int y = 0; y < 4; y++)
        for (int x = 0; x < 4; x++)
          rot[y][x] = values[x + 4 * y].toDouble();
      shot.Extrinsics.SetRot(rot);

      vcg::Camera<ScalarType> &cam = shot.Intrinsics;
      if (attr.contains("CameraType")) cam.cameraType = attr.namedItem("CameraType").nodeValue().toInt();
      cam.FocalMm = attr.namedItem("FocalMm").nodeValue().toDouble();
      cam.ViewportPx.X() = attr.namedItem("ViewportPx").nodeValue().section(' ', 0, 0).toInt();
      cam.ViewportPx.Y() = attr.namedItem("ViewportPx").nodeValue().section(' ', 1, 1).toInt();
      cam.CenterPx[0] = attr.namedItem("CenterPx").nodeValue().section(' ', 0, 0).toDouble();
      cam.CenterPx[1] = attr.namedItem("CenterPx").nodeValue().section(' ', 1, 1).toDouble();
      cam.PixelSizeMm[0] = attr.namedItem("PixelSizeMm").nodeValue().section(' ', 0, 0).toDouble();
      cam.PixelSizeMm[1] = attr.namedItem("PixelSizeMm").nodeValue().section(' ', 1, 1).toDouble();
      cam.k[0] = attr.namedItem("LensDistortion").nodeValue().section(' ', 0, 0).toDouble();
      cam.k[1] = attr.namedItem("LensDistortion").nodeValue().section(' ', 1, 1).toDouble();
    }
    // scale correction should no more exist !!!
    //    float scorr = attr.namedItem("ScaleCorr").nodeValue().toDouble();
    //    if(scorr != 0.0) {
    //      cam.PixelSizeMm[0] *= scorr;
    //      cam.PixelSizeMm[1] *= scorr;
    //    }
    return true;
  }
  return false;
}

/// TEXALIGN VERSION
template <class ShotType>
  bool ReadShotFromOLDXML( ShotType &shot, const QDomNode &node)
{
  typedef typename ShotType::ScalarType ScalarType;

  if(QString::compare(node.nodeName(),"CamParam")==0)
  {
    QDomNamedNodeMap attr = node.attributes();
    vcg::Point3<ScalarType> tra;
    tra[0] = attr.namedItem("SimTra").nodeValue().section(' ',0,0).toDouble();
    tra[1] = attr.namedItem("SimTra").nodeValue().section(' ',1,1).toDouble();
    tra[2] = attr.namedItem("SimTra").nodeValue().section(' ',2,2).toDouble();
    shot.Extrinsics.SetTra(-tra);

    vcg::Matrix44<ScalarType> rot;
    QStringList values =  attr.namedItem("SimRot").nodeValue().split(" ", QString::SkipEmptyParts);
    for(int y = 0; y < 4; y++)
      for(int x = 0; x < 4; x++)
        rot[y][x] = values[x + 4*y].toDouble();
    shot.Extrinsics.SetRot(rot);

    vcg::Camera<ScalarType> &cam = shot.Intrinsics;
    cam.FocalMm = attr.namedItem("Focal").nodeValue().toDouble();
    cam.ViewportPx.X() = attr.namedItem("Viewport").nodeValue().section(' ',0,0).toInt();
    cam.ViewportPx.Y() = attr.namedItem("Viewport").nodeValue().section(' ',1,1).toInt();
    cam.CenterPx[0] = attr.namedItem("Center").nodeValue().section(' ',0,0).toInt();
    cam.CenterPx[1] = attr.namedItem("Center").nodeValue().section(' ',1,1).toInt();
    cam.PixelSizeMm[0] = attr.namedItem("ScaleF").nodeValue().section(' ',0,0).toDouble();
    cam.PixelSizeMm[1] = attr.namedItem("ScaleF").nodeValue().section(' ',1,1).toDouble();
    cam.k[0] = attr.namedItem("LensDist").nodeValue().section(' ',0,0).toDouble();
    cam.k[1] = attr.namedItem("LensDist").nodeValue().section(' ',1,1).toDouble();

    // scale correction
    float scorr = attr.namedItem("ScaleCorr").nodeValue().toDouble();
    if(scorr != 0.0) {
      cam.PixelSizeMm[0] *= scorr;
      cam.PixelSizeMm[1] *= scorr;
    }
    return true;
  }
  return false;
}



template <class ShotType>
    QDomElement WriteShotToQDomNode(
        const ShotType &shot, /// the shot to be written node
        QDomDocument &doc) /// The XML node to be read
{
  typedef typename ShotType::ScalarType ScalarType;

  QDomElement shotElem = doc.createElement( "VCGCamera" );
  vcg::Point3<ScalarType> tra = -(shot.Extrinsics.Tra());
  QString str = QString("%1 %2 %3 1").arg(tra[0]).arg(tra[1]).arg(tra[2]);
  shotElem.setAttribute("TranslationVector", str);
  vcg::Matrix44<ScalarType> rot = shot.Extrinsics.Rot();
  str = QString("%1 %2 %3 %4 %5 %6 %7 %8 %9 %10 %11 %12 %13 %14 %15 %16 ")
    .arg(rot[0][0]).arg(rot[0][1]).arg(rot[0][2]).arg(rot[0][3])
    .arg(rot[1][0]).arg(rot[1][1]).arg(rot[1][2]).arg(rot[1][3])
    .arg(rot[2][0]).arg(rot[2][1]).arg(rot[2][2]).arg(rot[2][3])
    .arg(rot[3][0]).arg(rot[3][1]).arg(rot[3][2]).arg(rot[3][3]);
  shotElem.setAttribute( "RotationMatrix", str);

  const vcg::Camera<ScalarType> &cam = shot.Intrinsics;

  shotElem.setAttribute("CameraType", cam.cameraType);
  shotElem.setAttribute("BinaryData", 0);

  shotElem.setAttribute( "FocalMm", cam.FocalMm);

  str = QString("%1 %2").arg(cam.k[0]).arg(cam.k[1]);
  shotElem.setAttribute( "LensDistortion", str);

  str = QString("%1 %2").arg(cam.PixelSizeMm[0]).arg(cam.PixelSizeMm[1]);
  shotElem.setAttribute( "PixelSizeMm", str);

  str = QString("%1 %2").arg(cam.ViewportPx[0]).arg(cam.ViewportPx[1]);
  shotElem.setAttribute( "ViewportPx", str);

  str = QString("%1 %2").arg(cam.CenterPx[0]).arg(cam.CenterPx[1]);
  shotElem.setAttribute( "CenterPx", str);

  //scale correction should no more exist !!!
  //str = QString("%1").arg((double) 1);
  //shotElem.setAttribute( "ScaleCorr", str);
  return shotElem;
}


template <class ShotType>
QDomElement WriteShotToQDomNodeBinary(
  const ShotType &shot, /// the shot to be written node
  QDomDocument &doc) /// The XML node to be read
{
  typedef typename ShotType::ScalarType ScalarType;

  QDomElement shotElem = doc.createElement("VCGCamera");
  vcg::Point3<ScalarType> tra = -(shot.Extrinsics.Tra());

  QByteArray value = QByteArray::fromRawData((char *)tra.V(), sizeof(ScalarType) * 3).toBase64();
  shotElem.setAttribute("TranslationVector", QString(value));
  
  vcg::Matrix44<ScalarType> rot = shot.Extrinsics.Rot();
  value = QByteArray::fromRawData((char *)rot.V(), sizeof(ScalarType) * 16).toBase64();
  shotElem.setAttribute("RotationMatrix", QString(value));

  const vcg::Camera<ScalarType> &cam = shot.Intrinsics;

  shotElem.setAttribute("CameraType", cam.cameraType);

  shotElem.setAttribute("BinaryData", 1);

  value = QByteArray::fromRawData((char *)&cam.FocalMm, sizeof(ScalarType)).toBase64();
  shotElem.setAttribute("FocalMm", QString(value));

  value = QByteArray::fromRawData((char *)&cam.k, sizeof(ScalarType) * 2).toBase64();
  shotElem.setAttribute("LensDistortion", QString(value));

  value = QByteArray::fromRawData((char *)&cam.PixelSizeMm, sizeof(ScalarType) * 2).toBase64();
  shotElem.setAttribute("PixelSizeMm", QString(value));

  value = QByteArray::fromRawData((char *)&cam.ViewportPx, sizeof(int) * 2).toBase64();
  shotElem.setAttribute("ViewportPx", QString(value));

  //str = QString("%1 %2").arg(cam.CenterPx[0]).arg(cam.CenterPx[1]);
  value = QByteArray::fromRawData((char *)&cam.CenterPx, sizeof(ScalarType) * 2).toBase64();
  shotElem.setAttribute("CenterPx", QString(value));

  //scale correction should no more exist !!!
  //str = QString("%1").arg((double) 1);
  //shotElem.setAttribute( "ScaleCorr", str);
  return shotElem;
}


#endif // SHOT_QT_H
