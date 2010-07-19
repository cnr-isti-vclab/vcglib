#ifndef SHOT_QT_H
#define SHOT_QT_H


/** This function read a shot from a parsed XML node.

  */

template <class ShotType>
    bool ReadShotFromQDomNode(
        ShotType &shot, /// the shot that will contain the read node
        const QDomNode &node) /// The XML node to be read
{
  if(QString::compare(node.nodeName(),"CamParam")==0)
  {
    QDomNamedNodeMap attr = node.attributes();
    vcg::Point3d tra;
    tra[0] = attr.namedItem("SimTra").nodeValue().section(' ',0,0).toDouble();
    tra[1] = attr.namedItem("SimTra").nodeValue().section(' ',1,1).toDouble();
    tra[2] = attr.namedItem("SimTra").nodeValue().section(' ',2,2).toDouble();
    shot.Extrinsics.SetTra(-tra);

    vcg::Matrix44d rot;
    QStringList values =  attr.namedItem("SimRot").nodeValue().split(" ", QString::SkipEmptyParts);
    for(int y = 0; y < 4; y++)
      for(int x = 0; x < 4; x++)
        rot[y][x] = values[x + 4*y].toDouble();
    shot.Extrinsics.SetRot(rot);

    vcg::Camera<double> &cam = shot.Intrinsics;
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
#endif // SHOT_QT_H
