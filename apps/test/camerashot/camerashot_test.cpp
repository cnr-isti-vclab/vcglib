
// STD headers
#include <iostream>
#include <assert.h>

// VCG headers
#include <vcg/math/camera.h>
#include <vcg/math/shot.h>

double dist2(vcg::Point2d p1, vcg::Point2d p2)
{
	double d = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]));
	return d;
}

double dist3(vcg::Point3d p1, vcg::Point3d p2)
{
	double d = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));
	return d;
}

// TEST1 - PROJECT A 3D POINT IN WORLD COORDINATE ON THE IMAGE PLANE
///////////////////////////////////////////////////////////////////////////////
void test1(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	vcg::Point2d p1proj, p2proj;

	p1proj = shot1.Project(p1);
	p2proj = shot2.Project(p2);

	vcg::Point2d p1test(633.58101456110933, 400.0);
	vcg::Point2d p2test(289.02943695191425, 400.0);

	assert(dist2(p1proj,p1test) < 0.00000001);
	assert(dist2(p2proj,p2test) < 0.00000001);
}

// TEST 2 - PROJECTION AND UNPROJECTION
///////////////////////////////////////////////////////////////////////////////
void test2(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{	
	vcg::Point2d p1proj, p2proj;

	p1proj = shot1.Project(p1);
	p2proj = shot2.Project(p2);

	vcg::Point3d p1unproj, p2unproj;
	vcg::Point3d pcam1, pcam2;

	vcg::Matrix44d R1 = shot1.Extrinsics.Rot();
	vcg::Point4d pp(-10.0, -5.0, -70.0, 1.0);
	pp = R1 * pp;

	pcam1 = shot1.ConvertWorldToCameraCoordinates(p1);
	p1unproj = shot1.UnProject(p1proj, pcam1[2]);
	pcam2 = shot2.ConvertWorldToCameraCoordinates(p2);
	p2unproj = shot2.UnProject(p2proj, pcam2[2]);

	assert(dist3(p1, p1unproj) < 0.00000001);
	assert(dist3(p2, p2unproj) < 0.00000001);
}

// TEST 3 - CAMERA CONVERSION - CONVERT FOCAL IN PIXELS IN FOCAL IN MM
///////////////////////////////////////////////////////////////////////////////
void test3(vcg::Shotd shot, vcg::Point3d p1)
{
	vcg::Shotd shotpx = shot;

	// we assume focal is in pixels
	shotpx.Intrinsics.FocalMm = 1028.5949393128985805389837482;
	shotpx.Intrinsics.PixelSizeMm[0] = 1.0;
	shotpx.Intrinsics.PixelSizeMm[1] = 1.0;

	vcg::Point2d pproj;
	pproj = shotpx.Project(p1);

	vcg::Point2d p1proj;
	p1proj = shot.Project(p1);

	assert(dist2(pproj, p1proj) < 0.000000001);

	// CONVERSION - (ccd is assumed to be 35mm width)

	double ccd_width = 35.0; // conventionally assumed
	double ccd_height = (ccd_width * shotpx.Intrinsics.ViewportPx[1]) / shotpx.Intrinsics.ViewportPx[0];
	shotpx.Intrinsics.PixelSizeMm[0] = (ccd_width / shotpx.Intrinsics.ViewportPx[0]);
	shotpx.Intrinsics.PixelSizeMm[1] = (ccd_height / shotpx.Intrinsics.ViewportPx[1]);
	shotpx.Intrinsics.FocalMm = (ccd_width * shotpx.Intrinsics.FocalMm) / shotpx.Intrinsics.ViewportPx[0];  // NOW FOCAL IS IN MM

	p1proj = shotpx.Project(p1);

	assert(dist2(pproj, p1proj) < 0.000000001);
}

// TEST 4 - CAMERA MODIFICATION - CHANGE SCALE FACTOR OF THE WORLD
///////////////////////////////////////////////////////////////////////////////
void test4(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	//...TODO...
}

// TEST 5 - CAMERA MODIFICATION - ROTATION + TRANSLATION
///////////////////////////////////////////////////////////////////////////////
void test5(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	vcg::Matrix44d R = shot1.Extrinsics.Rot();
	vcg::Point3d tr = shot1.Extrinsics.Tra();

	vcg::Matrix44d T;
	T.SetIdentity();
	R.ElementAt(3,0) = tr[0];
	R.ElementAt(3,1) = tr[1];
	R.ElementAt(3,2) = tr[2];
	R.ElementAt(3,3) = 1.0;

	shot1.MultMatrix(T);
	shot1.MultMatrix(R.transpose());

	R = shot2.Extrinsics.Rot();
	tr = shot2.Extrinsics.Tra();

	T.ElementAt(3,0) = -tr[0];
	T.ElementAt(3,1) = -tr[1];
	T.ElementAt(3,2) = -tr[2];
	T.ElementAt(3,3) = 1.0;

	shot1.MultMatrix(R);
	shot1.MultMatrix(T);

	vcg::Point2d p1proj1, p2proj1, p1proj2, p2proj2;
	p1proj1 = shot1.Project(p1);
	p1proj2 = shot2.Project(p1);

	p2proj1 = shot1.Project(p2);
	p2proj2 = shot2.Project(p2);

	assert(dist2(p1proj1, p1proj2) < 0.00001);
	assert(dist2(p2proj1, p2proj2) < 0.00001);
}


int main() 
{
	vcg::Point3d p1(20.0, 25.0, 10.0);
	vcg::Point3d p2(-6.0, 40.0, 50.0);
	vcg::Shotd shot1;
	vcg::Shotd shot2;

	// Initialize camera 1 (C1)
	shot1.Intrinsics.cameraType = vcg::Camera<double>::PERSPECTIVE;
	shot1.Intrinsics.FocalMm = 30.0;
	shot1.Intrinsics.CenterPx[0] = 600.0; shot1.Intrinsics.CenterPx[1] = 400.0;
	shot1.Intrinsics.ViewportPx[0] = 1200; shot1.Intrinsics.ViewportPx[1] = 800;
	shot1.Intrinsics.PixelSizeMm[0] = 0.029166; shot1.Intrinsics.PixelSizeMm[1] = 0.029166;

	// no distorion is assumed (!)
	shot1.Intrinsics.DistorCenterPx[0] = shot1.Intrinsics.DistorCenterPx[1] = 0.0;
	shot1.Intrinsics.k[0] = 0.0; shot1.Intrinsics.k[1] = 0.0; shot1.Intrinsics.k[2] = 0.0; shot1.Intrinsics.k[3] = 0.0; 

	vcg::Matrix44d R1; // -10 degree around Y axis
	double deg2rad = 0.01745329251994329576923690768489;
	R1.ElementAt(0,0) = std::cos(-10.0*deg2rad);
	R1.ElementAt(0,1) = 0.0;
	R1.ElementAt(0,2) = std::sin(-10.0*deg2rad);
	R1.ElementAt(0,3) = 0.0;
	R1.ElementAt(1,0) = 0.0;
	R1.ElementAt(1,1) = 0.0;
	R1.ElementAt(1,2) = 0.0;
	R1.ElementAt(1,3) = 0.0;
	R1.ElementAt(2,0) = -std::sin(-10.0*deg2rad);
	R1.ElementAt(2,1) = 0.0;
	R1.ElementAt(2,2) = std::cos(-10.0*deg2rad);
	R1.ElementAt(2,3) = 0.0;
	R1.ElementAt(3,0) = 0.0;
	R1.ElementAt(3,1) = 0.0;
	R1.ElementAt(3,2) = 0.0;
	R1.ElementAt(3,3) = 1.0;

	vcg::Point3d T1(30.0, 30.0, 80.0);
	shot1.Extrinsics.SetTra(T1);
	shot1.Extrinsics.SetRot(R1);

	// Initialize camera 2 (C2)
	shot2.Intrinsics.cameraType = vcg::Camera<double>::PERSPECTIVE;
	shot2.Intrinsics.FocalMm = 30.0;
	shot2.Intrinsics.CenterPx[0] = 600.0; shot2.Intrinsics.CenterPx[1] = 400.0;
	shot2.Intrinsics.ViewportPx[0] = 1200; shot2.Intrinsics.ViewportPx[1] = 800;
	shot2.Intrinsics.PixelSizeMm[0] = 0.029166; shot2.Intrinsics.PixelSizeMm[1] = 0.029166;

	// no distortion is assumed (!)
	shot2.Intrinsics.DistorCenterPx[0] = shot2.Intrinsics.DistorCenterPx[1] = 0.0;
	shot2.Intrinsics.k[0] = 0.0; shot2.Intrinsics.k[1] = 0.0; shot2.Intrinsics.k[2] = 0.0; shot2.Intrinsics.k[3] = 0.0; 


	vcg::Matrix44d R2; // 18 degree around Y axis (+ 180 degree for the correct orientation of the camera)
	R2.ElementAt(0,0) = std::cos(-45.0*deg2rad);
	R2.ElementAt(0,1) = 0.0;
	R2.ElementAt(0,2) = std::sin(-45.0*deg2rad);
	R2.ElementAt(0,3) = 0.0;
	R2.ElementAt(1,0) = 0.0;
	R2.ElementAt(1,1) = 0.0;
	R2.ElementAt(1,2) = 0.0;
	R2.ElementAt(1,3) = 0.0;
	R2.ElementAt(2,0) = -std::sin(-45.0*deg2rad);
	R2.ElementAt(2,1) = 0.0;
	R2.ElementAt(2,2) = std::cos(-45.0*deg2rad);
	R2.ElementAt(2,3) = 0.0;
	R2.ElementAt(3,0) = 0.0;
	R2.ElementAt(3,1) = 0.0;
	R2.ElementAt(3,2) = 0.0;
	R2.ElementAt(3,3) = 1.0;

	vcg::Point3d T2(50.0, 30.0, 80.0);
	shot2.Extrinsics.SetTra(T2);
	shot2.Extrinsics.SetRot(R2);


	// TEST 1 - project a 3D point in World coordinates on the image plane
	test1(shot1, shot2, p1, p2);

	// TEST 2 - projection and unprojection
	test2(shot1, shot2, p1, p2);

	// TEST 3 - CAMERA CONVERSION - CONVERT FOCAL IN PIXELS IN FOCAL IN MM
	test3(shot1, p1);

	// TEST 4 - CAMERA MODIFICATION - CHANGE SCALE FACTOR OF THE WORLD
	test4(shot1, shot2, p1, p2);

	// TEST 5 - CAMERA MODIFICATION - ROTATION + TRANSLATION
	test5(shot1, shot2, p1, p2);

  return 0;
}
