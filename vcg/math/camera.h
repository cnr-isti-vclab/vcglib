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



#ifndef __VCGLIB_CAMERA
#define __VCGLIB_CAMERA

// VCG
#include <vcg/space/point3.h>
#include <vcg/space/point2.h>
#include <vcg/math/similarity.h>



namespace vcg{


template<class S>
class Camera
{
public:
    typedef S ScalarType;
    enum {
        PERSPECTIVE =	0,
        ORTHO       =	1,
        ISOMETRIC	=   2,
        CAVALIERI	=   3
    };

    Camera():
        FocalMm(0.f),
        ViewportPx(vcg::Point2<int>(0,0)),
        PixelSizeMm(vcg::Point2<S>(0.0,0.0)),
        CenterPx(vcg::Point2<S>(0.0,0.0)),
        DistorCenterPx(vcg::Point2<S>(0.0,0.0)),
        cameraType(0)
    {
        k[0] = k[1] = k[2] = k[3] = 0;
    }

    //------ camera intrinsics
    ScalarType	FocalMm;			/// Focal Distance: the distance between focal center and image plane. Expressed in mm
    Point2<int>	ViewportPx;			/// Dimension of the Image Plane (in pixels)
    Point2<S>	PixelSizeMm;		/// Dimension in mm of a single pixel
    Point2<S>	CenterPx;			/// Position of the projection of the focal center on the image plane. Expressed in pixels

    Point2<S>	DistorCenterPx;		/// Position of the radial distortion center on the image plane in pixels
    S			k[4];				/// 1st & 2nd order radial lens distortion coefficient (only the first 2 terms are used)
    //------------------------

    int cameraType;					/// Type of camera: PERSPECTIVE,ORTHO,ISOMETRIC,CAVALIERI

    void SetOrtho( S l,S r, S b, S t,  vcg::Point2<int> viewport)
    {
        cameraType = ORTHO;
        ViewportPx = viewport;

        PixelSizeMm[0] = (r-l) / (S)ViewportPx[0];
        PixelSizeMm[1] = (t-b) / (S)ViewportPx[1];

        CenterPx[0] = -l/(r-l) * (S)ViewportPx[0];
        CenterPx[1] = -b/(t-b) * (S)ViewportPx[1];

    };


    bool IsOrtho() const
    {
        return ( cameraType == ORTHO );
    }

    //--- Set-up methods

    /// set the camera specifying the perspecive view
    inline void SetPerspective(S AngleDeg, S AspectRatio, S Focal, vcg::Point2<int> Viewport);

    /// set the camera specifying the cavalieri view
    inline void SetCavalieri(S sx, S dx, S bt, S tp, S Focal, vcg::Point2<int> Viewport);

    /// set the camera specifying the isometric view
    inline void SetIsometric(S sx, S dx, S bt, S tp, S Focal, vcg::Point2<int> Viewport);

    /// set the camera specifying the frustum view
    inline void SetFrustum(S dx, S sx, S bt, S tp, S Focal, vcg::Point2<int> Viewport);
    //------------------

    /// returns the projection matrix
    vcg::Matrix44<S> GetMatrix(S nearVal, S farVal);

    /// returns the frustum
    inline void GetFrustum(S & sx, S & dx, S & bt, S & tp, S & nr);

    //--- Space transformation methods

    /// project a point from 3d CAMERA space to the camera local plane
    inline vcg::Point2<S> Project(const vcg::Point3<S> & p) const;

    /// unproject a point from the camera local plane (plus depth) to 3d CAMERA space
    inline vcg::Point3<S> UnProject(const vcg::Point2<S> & p, const S & d) const;

    /// transforms local plane coords to vieport (pixel) coords
    inline vcg::Point2<S> LocalToViewportPx(const vcg::Point2<S> & p) const;

    /// transforms vieport (pixel) coords to local plane coords
    inline vcg::Point2<S> ViewportPxToLocal(const vcg::Point2<S> & p) const;

    /// transforms vieport (pixel) coords to [-1 1] coords
    inline vcg::Point2<S> ViewportPxTo_neg1_1(const vcg::Point2<S> & p) const;

    /// transforms [-1 1] coords to vieport (pixel) coords MICHELE IO
    inline vcg::Point2<S> Neg1_1ToViewportPx(const vcg::Point2<S> & p) const;

    /// transforms local plane coords to [0 1] coords
    inline vcg::Point2<S> LocalTo_0_1(const vcg::Point2<S> & p) const;

    /// transforms local plane coords to [-1 1] coords
    inline vcg::Point2<S> LocalTo_neg1_1(const vcg::Point2<S> & p) const;

    /// transforms an undistorted 2D camera plane point in a distorted 2D camera plane point
    vcg::Point2<S> UndistortedToDistorted(vcg::Point2<S> u) const;

    /// transforms a distorted 2D camera plane point in an undistorted 2D camera plane point
    vcg::Point2<S> DistortedToUndistorted(vcg::Point2<S> d) const;
    //--------------------------------

	/// import from a different type
	template <class S1> 
	Camera<S> Import(const Camera<S1> & cam){
		this->FocalMm = cam.FocalMm;			
		this->ViewportPx[0] = cam.ViewportPx[0];			
		this->ViewportPx[1] = cam.ViewportPx[1];			
		this->PixelSizeMm.Import(cam.PixelSizeMm);		
		this->CenterPx.Import(cam.CenterPx);			
		this->DistorCenterPx.Import(cam.DistorCenterPx);
		this->k[0] = cam.k[0];				
		this->k[1] = cam.k[1];				
		this->k[2] = cam.k[2];				
		this->k[3] = cam.k[3];		
		this->cameraType = cam.cameraType;
		return *this;
	};
};



/// project a point from 3d CAMERA space to the camera's plane]
template<class S>
vcg::Point2<S> Camera<S>::Project(const vcg::Point3<S> & p) const
{
    vcg::Point2<S> q =  Point2<S>(p[0],p[1]);
    vcg::Point2<S> d =  Point2<S>(p[0],p[1]);

    if(!IsOrtho())
    {
        q[0] *= FocalMm/p.Z();
        q[1] *= FocalMm/p.Z();

        if(k[0]!=0)
        {
            vcg::Point2<S> d;
            d=UndistortedToDistorted(q);
            q=d;
        }

    }


    return q;
}

template <class S> vcg::Matrix44<S> Camera<S>::GetMatrix(S nearVal, S farVal) {
    S left, right, bottom, top, nr;
    GetFrustum(left, right, bottom,top, nr);

    if(cameraType == PERSPECTIVE) {
      S ratio = nearVal/nr;
      left *= ratio;
      right *= ratio;
      bottom *= ratio;
      top *= ratio;
    }
    vcg::Matrix44<S> m;
    m[0][0] = 2.0*nearVal/(right - left);
    m[0][1] = 0;
    m[0][2] = (right + left)/(right - left);
    m[0][3] = 0;

    m[1][0] = 0;
    m[1][1] = 2*nearVal/(top - bottom);
    m[1][2] = (top + bottom)/(top - bottom);
    m[1][3] = 0;

    m[2][0] = 0;
    m[2][1] = 0;
    m[2][2] = -(farVal + nearVal)/(farVal - nearVal);
    m[2][3] =  - 2*farVal*nearVal/(farVal - nearVal);

    m[3][0] = 0;
    m[3][1] = 0;
    m[3][2] = -1;
    m[3][3] = 0;

    return m;
}

/// unproject a point from the camera 2d plane [-1,-1]x[1,1] (plus depth) to 3d CAMERA space
template<class S>
vcg::Point3<S> Camera<S>::UnProject(const vcg::Point2<S> & p, const S & d) const
{
    vcg::Point3<S> np = Point3<S>(p[0], p[1], d);

    if(!IsOrtho())
    {
        if(k[0]!=0)
        {
            vcg::Point2<S> d = Point2<S>(p[0], p[1]);
            vcg::Point2<S> u=DistortedToUndistorted(d);
            np[0]=u[0];
            np[1]=u[1];
        }

        np[0] /= FocalMm/d;
        np[1] /= FocalMm/d;
    }

    return np;
}

/// transforms local plane coords to vieport (pixel) coords
template<class S>
vcg::Point2<S> Camera<S>::LocalToViewportPx(const vcg::Point2<S> & p) const
{
    vcg::Point2<S> np;

    np[0] = (p[0] / PixelSizeMm.X()) + CenterPx.X();
    np[1] = (p[1] / PixelSizeMm.Y()) + CenterPx.Y();

    return np;
}

/// transforms vieport (pixel) coords to local plane coords
template<class S>
vcg::Point2<S> Camera<S>::ViewportPxToLocal(const vcg::Point2<S> & p) const
{
    vcg::Point2<S>  ps;
    ps[0] = (p[0]-CenterPx.X()) * PixelSizeMm.X();
    ps[1] = (p[1]-CenterPx.Y()) * PixelSizeMm.Y();
    return ps;
}

/// transforms vieport (pixel) coords to [-1 1] coords
template<class S>
vcg::Point2<S> Camera<S>::ViewportPxTo_neg1_1(const vcg::Point2<S> & p) const
{
    vcg::Point2<S>  ps;
    ps[0] = 2.0f * ((p[0]-CenterPx.X()) * PixelSizeMm.X())/ ( PixelSizeMm.X() * (S)ViewportPx[0] );
    ps[1] = 2.0f * ((p[1]-CenterPx.Y()) * PixelSizeMm.Y())/ ( PixelSizeMm.Y() * (S)ViewportPx[1] );
    return ps;
}

/// transforms [-1 1] coords to vieport (pixel) coords MICHELE IO
template<class S>
vcg::Point2<S> Camera<S>::Neg1_1ToViewportPx(const vcg::Point2<S> & p) const
{
    vcg::Point2<S>  ps;
    ps[0] = ((PixelSizeMm.X() * (S)ViewportPx[0] *p[0])/(2.0f * PixelSizeMm.X()))+CenterPx.X();
    ps[1] = ((PixelSizeMm.Y() * (S)ViewportPx[1] *p[1])/(2.0f * PixelSizeMm.Y()))+CenterPx.Y();
    return ps;
}

/// transforms local plane coords to [0-1] coords
template<class S>
vcg::Point2<S> Camera<S>::LocalTo_0_1(const vcg::Point2<S> & p) const
{
    vcg::Point2<S>  ps;
    ps[0] = ( p[0]/PixelSizeMm.X() + CenterPx.X() ) / (S)ViewportPx[0];
    ps[1] = ( p[1]/PixelSizeMm.Y() + CenterPx.Y() ) / (S)ViewportPx[1];
    return ps;
}

/// transforms local plane coords to [-1 1] coords
template<class S>
vcg::Point2<S> Camera<S>::LocalTo_neg1_1(const vcg::Point2<S> & p) const
{
    vcg::Point2<S>  ps;
    ps[0] = 2.0f * p[0] / ( PixelSizeMm.X() * (S)ViewportPx[0] );
    ps[1] = 2.0f * p[1] / ( PixelSizeMm.Y() * (S)ViewportPx[1] );
    return ps;
}

/// transforms an undistorted 2D camera plane point in a distorted 2D camera plane point
template<class Scalar>
vcg::Point2<Scalar> Camera<Scalar>::UndistortedToDistorted(vcg::Point2<Scalar>  u) const
{
    vcg::Point2<Scalar> dis;
    vcg::Point2<Scalar> dc=ViewportPxTo_neg1_1(DistorCenterPx);
    const Scalar SQRT3 = Scalar(1.732050807568877293527446341505872366943);
    const Scalar CBRT = Scalar(0.33333333333333333333333);
    Scalar Ru,Rd,lambda,c,d,Q,R,D,S,T,sinT,cosT;

    if(((u[0]-dc[0])==0 && (u[1]-dc[1])==0) || k[0] == 0)
    {
        dis[0] = u[0];
        dis[1] = u[1];
        return dis;
    }

    Ru = hypot ((u[0]-dc[0]), (u[1]-dc[1]));	/* SQRT(Xu*Xu+Yu*Yu) */
    c = 1 / k[0];
    d = -c * Ru;

    Q = c / 3;
    R = -d / 2;
    if (R<0)
        D = pow(Q,3) + sqrt(-R);
    else
        D = pow(Q,3) + sqrt(R);

    if (D >= 0)		/* one real root */
    {
        D = sqrt(D);
        S = pow((R + D),CBRT);
        if (R>=D)
            T = pow((R - D),CBRT);
        else
            T = - pow(abs((int)(R - D)),CBRT); //MODIFICATO DA ME
        Rd = S + T;

        if (Rd < 0)
            Rd = sqrt(-1 / (3 * k[0]));
    }
    else			/* three real roots */
    {
        D = sqrt(-D);
        S = pow((Scalar)(hypot (R, D)),(Scalar)CBRT);
        T = atan2 (D, R) / 3;
        //SinCos(T, sinT, cosT);
        sinT=sin(T);
        cosT=cos(T);

        /* the larger positive root is    2*S*cos(T)                   */
        /* the smaller positive root is   -S*cos(T) + SQRT(3)*S*sin(T) */
        /* the negative root is           -S*cos(T) - SQRT(3)*S*sin(T) */
        Rd = -S * cosT + SQRT3 * S * sinT;	/* use the smaller positive root */
    }

    lambda = Rd / Ru;

    dis[0] = u[0] * lambda;
    dis[1] = u[1] * lambda;

    return dis;
}

/// transforms a distorted 2D camera plane point in an undistorted 2D camera plane point
template<class S>
vcg::Point2<S> Camera<S>::DistortedToUndistorted(vcg::Point2<S> d) const
{
    vcg::Point2<S> u;
    vcg::Point2<S> dc=ViewportPxTo_neg1_1(DistorCenterPx);
    S r=sqrt(pow((d[0]-dc[0]),2)+pow((d[1]-dc[1]),2));
    u[0]=d[0]*(1-k[0]*r*r);
    u[1]=d[1]*(1-k[0]*r*r);

    return u;

}


//--- basic camera setup (GL-like)

/// set the camera specifying the perspective view
template<class S>
void Camera<S>::SetPerspective(	S AngleDeg, S AspectRatio, S Focal, vcg::Point2<int> Viewport)
{
    cameraType = PERSPECTIVE;
    S halfsize[2];

    halfsize[1] = tan(math::ToRad(AngleDeg/2.0f)) * Focal;
    halfsize[0] = halfsize[1]*AspectRatio;

    SetFrustum(-halfsize[0],halfsize[0],-halfsize[1],halfsize[1],Focal,Viewport);
}

/// set the camera specifying the frustum view
template<class S>
void Camera<S>::SetFrustum(	S sx, S dx, S bt, S tp, S Focal, vcg::Point2<int> Viewport)
{
    S vp[2];
    vp[0] = dx-sx;
    vp[1] = tp-bt;

    ViewportPx[0] = Viewport[0];
    if(vp[1] != -1)
        ViewportPx[1] = Viewport[1];			// the user specified the viewport
    else
        ViewportPx[1] = ViewportPx[0];	// default viewport

    PixelSizeMm[0] = vp[0] / (S)Viewport[0];
    PixelSizeMm[1] = vp[1] / (S)Viewport[1];

    CenterPx[0] = -sx/vp[0] * (S)Viewport[0];
    CenterPx[1] = -bt/vp[1] * (S)Viewport[1];

    FocalMm =Focal;
}

//--- special cameras setup

/// set the camera specifying the cavalieri view
template<class S>
void Camera<S>::SetCavalieri(S sx, S dx, S bt, S tp, S Focal, vcg::Point2<int> Viewport)
{
    cameraType = CAVALIERI;
    SetFrustum(sx, dx, bt, tp, Focal,Viewport);
}

/// set the camera specifying the isometric view
template<class S>
void Camera<S>::SetIsometric(S sx, S dx, S bt, S tp, S Focal, vcg::Point2<int> Viewport)
{
    cameraType = ISOMETRIC;
    SetFrustum(sx, dx, bt, tp, Focal,Viewport);
}

/// returns the frustum
template<class S>
void Camera<S>:: GetFrustum( S & sx, S & dx, S & bt, S & tp, S & nr)
{
    dx = CenterPx.X()* PixelSizeMm.X();			//scaled center
    sx = -( (S)ViewportPx.X() - CenterPx.X() ) * PixelSizeMm.X();

    bt = -CenterPx.Y()* PixelSizeMm.Y();
    tp = ( (S)ViewportPx.Y() - CenterPx.Y() ) * PixelSizeMm.Y();

    nr = FocalMm;
}

}


#endif




