/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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
Revision 1.12  2005/07/11 13:12:35  cignoni
small gcc-related compiling issues (typenames,ending cr, initialization order)

Revision 1.11  2005/01/22 11:20:20  ponchio
<...Point3.h> -> <...point3.h>

Revision 1.10  2005/01/05 13:26:15  ganovelli
corretto cambiamento di sistema di rif.

Revision 1.9  2004/12/15 18:45:50  tommyfranken
*** empty log message ***

Revision 1.4  2004/10/07 14:41:31  fasano
Little fix on ViewPoint() method

Revision 1.3  2004/10/07 14:24:53  ganovelli
added LookAt,LookToward

Revision 1.2  2004/10/05 19:04:25  ganovelli
version 5-10-2004 in progress

Revision 1.1  2004/09/15 22:58:05  ganovelli
re-creation

Revision 1.2  2004/09/06 21:41:30  ganovelli
*** empty log message ***

Revision 1.1  2004/09/03 13:01:51  ganovelli
creation

****************************************************************************/


#ifndef __VCGLIB_SHOT
#define __VCGLIB_SHOT

// #include <vector>
// #include <vcg/Matrix44.h>
// #include <vcg/Box3.h>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>
#include <vcg/math/similarity.h>
#include <vcg/math/camera.h>

namespace vcg{

template <class S>
class Shot {
public:
	typedef Camera<S> CameraType;
	typedef S ScalarType;

protected:	

	enum {
		NOTVALID_BIT = 0x1,
		CREATED_EMPTY= 0x2
	};
	char flags;
	Camera<S> & camera;																	// the camera that shot
	vcg::Similarity<S,vcg::Matrix44<S> > similarity;		// the position from where it did it

public:

	Shot( Camera<S> & c):camera(c){similarity.SetIdentity();}
	Shot():camera(*new vcg::Camera<S>()){similarity.SetIdentity();flags|=CREATED_EMPTY;}
	~Shot(){if(flags&CREATED_EMPTY) delete &camera;}

	bool IsValid(){ return (flags& NOTVALID_BIT)==0;}
	void SetValid(bool v){ if(!v) flags|=NOTVALID_BIT; else flags&=~NOTVALID_BIT;}

	Camera<S> & Camera(){return camera;};

	/// take the i-th axis of the coordinate system of the camera
	vcg::Point3<S> Axis(const int & i)const;

	/// take the viewpoint
	 const vcg::Point3<S> ViewPoint()const;

	/// set the viewpoint
	void SetViewPoint(const vcg::Point3<S> & viewpoint);

	/// look at
	void LookAt(const vcg::Point3<S> & z_dir,const vcg::Point3<S> & up);

	/// look at (same as opengl)
	void LookAt(const S & eye_x,const S & eye_y,const S & eye_z,const S & at_x,const S & at_y,const S & at_z,
							const S & up_x,const S & up_y,const S & up_z);

	/// look towards
	void LookTowards(const vcg::Point3<S> & z_dir,const vcg::Point3<S> & up);

	/// convert a 3d point in camera coordinates
	vcg::Point3<S>  ConvertToCameraCoordinates(const vcg::Point3<S> & p) const;

	/// convert a 3d point in camera coordinates
	vcg::Point3<S>  ConvertToWorldCoordinates(const vcg::Point3<S> & p) const;

	/// project onto the camera plane
	vcg::Point2<S> Project(const vcg::Point3<S> & p) const;

	vcg::Point3<S> UnProject(const vcg::Point2<S> & p) const;

	/// take the distance from the point p and the plane parallel to the camera plane and passing through the view
	/// point. The would be z depth 
	S Depth(const vcg::Point3<S> & p)const;

}; // end class definition

template <class S>
const vcg::Point3<S> Shot<S>::ViewPoint() const {
	//Matrix44<S> m = similarity.Matrix();
	//return Point3<S>(m[0][3],m[1][3],m[2][3]);
	return -similarity.tra;
}

template <class S>
	vcg::Point3<S>  Shot<S>::Axis(const int & i) const {	
			vcg::Matrix44<S> m; 
			similarity.rot.ToMatrix(m); 
			vcg::Point3<S> aa = m.Row3(i);
			return aa;
	}

template <class S>
void Shot<S>::SetViewPoint(const vcg::Point3<S> & viewpoint){
	similarity.tra = -viewpoint;
}

template <class S>
void Shot<S>::LookAt(const vcg::Point3<S> & z_dir,const vcg::Point3<S> & up){
	  LookTowards(z_dir-ViewPoint(),up);
}

template <class S>
void Shot<S>::LookAt(const S & eye_x,const S & eye_y,const S & eye_z,const S & at_x,const S & at_y,const S & at_z,
										 const S & up_x,const S & up_y,const S & up_z){
	SetViewPoint(Point3<S>(eye_x,eye_y,eye_z));
	LookAt(Point3<S>(at_x,at_y,at_z),Point3<S>(up_x,up_y,up_z));
}


template <class S>
void Shot<S>::LookTowards(const vcg::Point3<S> & z_dir,const vcg::Point3<S> & up){
	vcg::Point3<S> x_dir = up ^-z_dir ;
	vcg::Point3<S> y_dir = -z_dir ^x_dir ;
	
	Matrix44<S> m;
	m.SetIdentity();
	*(vcg::Point3<S> *)&m[0][0] = x_dir/x_dir.Norm();
	*(vcg::Point3<S> *)&m[1][0] = y_dir/y_dir.Norm();
	*(vcg::Point3<S> *)&m[2][0] = -z_dir/z_dir.Norm();

	similarity.rot.FromMatrix(m);
}

template <class S>
vcg::Point3<S> Shot<S>::ConvertToCameraCoordinates(const vcg::Point3<S> & p) const{
	Matrix44<S> rotM;
	similarity.rot.ToMatrix(rotM);
	vcg::Point3<S> cp = rotM * (p+similarity.tra);
	// note: the camera reference system is right handed
	cp[2]=-cp[2];
	return cp;
	}
template <class S>
vcg::Point3<S> Shot<S>::ConvertToWorldCoordinates(const vcg::Point3<S> & p) const{
	Matrix44<S> rotM;
	vcg::Point3<S> cp = p;
	// note: the World reference system is left handed
	cp[2]=-cp[2];
	similarity.rot.ToMatrix(rotM);
	cp = Inverse(rotM) * cp-similarity.tra;
	return cp;
}
template <class S>
vcg::Point2<S> Shot<S>::Project(const vcg::Point3<S> & p) const{
		return camera.Project(ConvertToCameraCoordinates(p));
	}
template <class S>
vcg::Point3<S> Shot<S>::UnProject(const vcg::Point2<S> & p) const{
	vcg::Point3<S> q = camera.UnProject(p);
	return ConvertToWorldCoordinates(q);
}
template <class S>
S Shot<S>::Depth(const vcg::Point3<S> & p)const {
	return ConvertToCameraCoordinates(p).Z();
}


class Shotf: public Shot<float>{};
class Shotd: public Shot<double>{};

};

#endif




