
// STD headers
#include <iostream>
#include <assert.h>

// VCG headers
#include <vcg/math/camera.h>
#include <vcg/math/shot.h>

static double precision = 0.000000001;  // 1e-9

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

double checkIdentity(vcg::Matrix44d M)
{
	double d;

	d = (M.ElementAt(0,0) - 1.0) * (M.ElementAt(0,0) - 1.0) + M.ElementAt(0,1) * M.ElementAt(0,1) + M.ElementAt(0,2) * M.ElementAt(0,2) + M.ElementAt(0,3) * M.ElementAt(0,3);
	d += M.ElementAt(1,0) * M.ElementAt(1,0) + (M.ElementAt(1,1) - 1.0) * (M.ElementAt(1,1) - 1.0) + M.ElementAt(1,2) * M.ElementAt(1,2) + M.ElementAt(1,3) * M.ElementAt(1,3);
	d += M.ElementAt(2,0) * M.ElementAt(2,0) + M.ElementAt(2,1) * M.ElementAt(2,1) + (M.ElementAt(2,2) - 1.0) * (M.ElementAt(2,2) - 1.0) + M.ElementAt(2,3) * M.ElementAt(2,3);
	d += M.ElementAt(3,0) * M.ElementAt(3,0) + M.ElementAt(3,1) * M.ElementAt(3,1) + M.ElementAt(3,2) * M.ElementAt(3,2) + (M.ElementAt(3,3) - 1.0) * (M.ElementAt(3,3) - 1.0);

	return d;
}

// TEST1 - PROJECT A 3D POINT IN WORLD COORDINATE ON THE IMAGE PLANE
///////////////////////////////////////////////////////////////////////////////
bool test1(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	vcg::Point2d p1proj, p2proj;

	p1proj = shot1.Project(p1);
	p2proj = shot2.Project(p2);

	vcg::Point2d p1test(633.58101456110933, 327.22860336234237);
	vcg::Point2d p2test(289.02943695191425, 315.42715619973069);

	if (dist2(p1proj,p1test) > precision)
		return false;

	if (dist2(p2proj,p2test) > precision)
		return false;

	return true;
}

// TEST 2 - PROJECTION AND UNPROJECTION
///////////////////////////////////////////////////////////////////////////////
bool test2(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{	
	vcg::Point2d p1proj, p2proj;

	p1proj = shot1.Project(p1);
	p2proj = shot2.Project(p2);

	vcg::Point3d p1unproj, p2unproj;
	vcg::Point3d pcam1, pcam2;

	pcam1 = shot1.ConvertWorldToCameraCoordinates(p1);
	p1unproj = shot1.UnProject(p1proj, pcam1[2]);
	pcam2 = shot2.ConvertWorldToCameraCoordinates(p2);
	p2unproj = shot2.UnProject(p2proj, pcam2[2]);

	if (dist3(p1, p1unproj) > precision)
		return false;

	if (dist3(p2, p2unproj) > precision)
		return false;

	return true;
}

// TEST 3 - DEPTH COMPUTATION
///////////////////////////////////////////////////////////////////////////////
bool test3(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	vcg::Point2d p1proj, p2proj;

	p1proj = shot1.Project(p1);
	p2proj = shot2.Project(p2);

	vcg::Point3d p1unproj, p2unproj;

	double depth1 = shot1.Depth(p1);
	double depth2 = shot2.Depth(p2);

	p1unproj = shot1.UnProject(p1proj, depth1);
	p2unproj = shot2.UnProject(p2proj, depth2);

	if (dist3(p1, p1unproj) > precision)
		return false;

	if (dist3(p2, p2unproj) > precision)
		return false;

	return true;
}

// TEST 4 - CAMERA CONVERSION - CONVERT FOCAL IN PIXELS IN FOCAL IN MM
///////////////////////////////////////////////////////////////////////////////
bool test4(vcg::Shotd shot, vcg::Point3d p1)
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

	if(dist2(pproj, p1proj) > precision)
		return false;

	// CONVERSION - (ccd is assumed to be 35mm width)
	shot.ConvertFocalToMM();

	p1proj = shotpx.Project(p1);

	if (dist2(pproj, p1proj) > precision)
		return false;

	return true;
}

// TEST 5 - CAMERA-SHOT MODIFICATION - CHANGE SCALE FACTOR OF THE WORLD
///////////////////////////////////////////////////////////////////////////////
bool test5(vcg::Shotd shot, vcg::Point3d p1, vcg::Point3d p2)
{
	vcg::Point2d p1projPX, p2projPX;
	p1projPX = shot.Project(p1);
	p2projPX = shot.Project(p2);

	vcg::Point2d p1proj, p2proj;
	p1proj = shot.Intrinsics.ViewportPxToLocal(p1projPX);
	p2proj = shot.Intrinsics.ViewportPxToLocal(p2projPX);

	vcg::Point2d diff;
	double distance_before_world_scaling;
	diff = (p2proj-p1proj);
  distance_before_world_scaling = diff.Norm();

	// WORLD SCALING
	double scalefactor = 20.0;

	// adjust World 3D points
	p1 *= scalefactor;
	p2 *= scalefactor;

	shot.RescalingWorld(scalefactor);

	p1projPX = shot.Project(p1);
	p2projPX = shot.Project(p2);

	p1proj = shot.Intrinsics.ViewportPxToLocal(p1projPX);
	p2proj = shot.Intrinsics.ViewportPxToLocal(p2projPX);

	double distance_after_world_scaling;
	diff = (p2proj - p1proj);
	distance_after_world_scaling = diff.Norm();

	if (std::fabs(distance_before_world_scaling - (distance_after_world_scaling / scalefactor)) > precision)
		return false;

	return true;
}

// TEST 6 - WORLD-TO-EXTRINSICS AND EXTRINSICS-TO-WORLD TRANSFORMATIONS
bool test6(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	vcg::Matrix44d WtoE1 = shot1.GetWorldToExtrinsicsMatrix();
	vcg::Matrix44d WtoE2 = shot2.GetWorldToExtrinsicsMatrix();
	vcg::Matrix44d EtoW1 = shot1.GetExtrinsicsToWorldMatrix();
	vcg::Matrix44d EtoW2 = shot2.GetExtrinsicsToWorldMatrix();

	vcg::Matrix44d I1 = WtoE1 * EtoW1;
	vcg::Matrix44d I2 = WtoE2 * EtoW2;
	vcg::Matrix44d I3 = EtoW1 * WtoE1;
	vcg::Matrix44d I4 = EtoW2 * WtoE2;

	if (checkIdentity(I1) > precision)
		return false;

	if (checkIdentity(I2) > precision)
		return false;

	if (checkIdentity(I3) > precision)
		return false;

	if (checkIdentity(I4) > precision)
		return false;

	vcg::Point3d axisX(1.0, 0.0, 0.0);
	vcg::Point3d axisY(0.0, 1.0, 0.0);
	vcg::Point3d axisZ(0.0, 0.0, 1.0);

	vcg::Point3d vx = EtoW1 * axisX;
	vcg::Point3d vy = EtoW1 * axisY;
	vcg::Point3d vz = EtoW1 * axisZ;

	if (dist3(vx, shot1.Extrinsics.Tra() + shot1.Extrinsics.Rot().GetRow3(0)) > precision)
		return false;
	
	if (dist3(vy, shot1.Extrinsics.Tra() + shot1.Extrinsics.Rot().GetRow3(1)) > precision)
		return false;
	
	if (dist3(vz, shot1.Extrinsics.Tra() + shot1.Extrinsics.Rot().GetRow3(2)) > precision)
		return false;

	return true;
}


// TEST 7 - SHOT MODIFICATION - ROTATION + TRANSLATION
///////////////////////////////////////////////////////////////////////////////
bool test7(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	vcg::Matrix44d R;
	R.SetZero();
	R.ElementAt(0,2) = 1.0;
	R.ElementAt(1,1) = 1.0;
	R.ElementAt(2,0) = -1.0;
	R.ElementAt(3,3) = 1.0;

	vcg::Point2d p1proj = shot1.Project(p1);

	vcg::Point3d prot = R * p1;
	shot1.ApplyRigidTransformation(R);
	vcg::Point2d protproj = shot1.Project(prot);

	return true;
}

// TEST 8 - SHOT MODIFICATION - ROTATION + TRANSLATION
///////////////////////////////////////////////////////////////////////////////
bool test8(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	// put shot1 reference frame into the origin of the World coordinates system
	vcg::Matrix44d M = shot1.GetWorldToExtrinsicsMatrix();
	shot1.ApplyRigidTransformation(M);

	// then, put in the shot2 reference frame
	M = shot2.GetExtrinsicsToWorldMatrix();
	shot1.ApplyRigidTransformation(M);

	// test..
	vcg::Point2d p1proj1, p2proj1, p1proj2, p2proj2;
	p1proj1 = shot1.Project(p1);
	p1proj2 = shot2.Project(p1);

	p2proj1 = shot1.Project(p2);
	p2proj2 = shot2.Project(p2);

	if (dist2(p1proj1, p1proj2) > precision)
		return false;

	if (dist2(p2proj1, p2proj2) > precision)
		return false;

	return true;
}

// TEST 9 - SHOT MODIFICATION - ROTATION + TRANSLATION
///////////////////////////////////////////////////////////////////////////////
bool test9(vcg::Shotd shot1, vcg::Shotd shot2, vcg::Point3d p1, vcg::Point3d p2)
{
	vcg::Matrix44d M1 = shot1.GetExtrinsicsToWorldMatrix();
	vcg::Matrix44d M2 = shot2.GetWorldToExtrinsicsMatrix();
	vcg::Matrix44d M;
	M = M2 * M1;  // roto-translation that maps the frame of Shot1 in the frame of Shot2

	// apply it..
	shot2.ApplyRigidTransformation(M);

	// and test it..
	vcg::Point2d p1proj1, p2proj1, p1proj2, p2proj2;
	p1proj1 = shot1.Project(p1);
	p1proj2 = shot2.Project(p1);

	p2proj1 = shot1.Project(p2);
	p2proj2 = shot2.Project(p2);

	if (dist2(p1proj1, p1proj2) > precision)
		return false;

	if (dist2(p2proj1, p2proj2) > precision)
		return false;

	return true;
}

int main() 
{
	vcg::Point3d p1(20.0, 25.0, 10.0);
	vcg::Point3d p2(-6.0, 25.0, 50.0);
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
	R1.ElementAt(1,1) = 1.0;
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
	R2.ElementAt(1,1) = 1.0;
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
	if (test1(shot1, shot2, p1, p2))
	{
		std::cout << "TEST 1 (projection) - PASSED(!)" << std::endl;
	}
	else 
		std::cout << "TEST 1 (projection) - FAILED(!)" << std::endl;

	// TEST 2 - projection and unprojection
	if (test2(shot1, shot2, p1, p2))
	{
		std::cout << "TEST 2 (unprojection) - PASSED(!)" << std::endl;
	}
	else
	{
		std::cout << "TEST 2 (unprojection) - FAILED(!)" << std::endl;
	}

	// TEST 3 - DEPTH COMPUTATION
	if (test3(shot1, shot2, p1, p2))
	{
		std::cout << "TEST 3 (depth computation) - PASSED(!)" << std::endl;
	}
	else
	{
		std::cout << "TEST 3 (depth computation) - FAILED(!)" << std::endl;
	}


	// TEST 4 - CAMERA CONVERSION - CONVERT FOCAL IN PIXELS IN FOCAL IN MM
	if (test4(shot1, p1))
	{
		std::cout << "TEST 4 (focal in px to focal in mm) - PASSED(!)" << std::endl;
	}
	else
	{
		std::cout << "TEST 4 (focal in px to focal in mm) - FAILED(!)" << std::endl;
	}

	// TEST 5 - CAMERA-SHOT MODIFICATION - CHANGE SCALE FACTOR OF THE WORLD
	if (test5(shot1, p1, p2))
	{
		std::cout << "TEST 5 (scaling the World) - PASSED(!)" << std::endl;
	}
	else
	{
		std::cout << "TEST 5 (scaling the World) - FAILED(!)" << std::endl;
	}

	// TEST 6 - WORLD-TO-EXTRINSICS AND EXTRINSICS-TO-WORLD TRANSFORMATIONS
	if (test6(shot1, shot2, p1, p2))
	{
		std::cout << "TEST 6 (World-to-Extrinsics and Extrinsics-to-World) - PASSED(!)" << std::endl;
	}
	else
	{
		std::cout << "TEST 6 (World-to-Extrinsics and Extrinsics-to-World) - FAILE(!)" << std::endl;
	}

	// TEST 7 - SHOT MODIFICATION - ROTO-TRANSLATION OF THE SHOT COORDINATES SYSTEM
	if (test7(shot1, shot2, p1, p2))
	{
		std::cout << "TEST 7 (roto-translation of the Shot coordinates system) - PASSED(!)" << std::endl;
	}
	else
	{
		std::cout << "TEST 7 (roto-translation of the Shot coordinates system) - FAILED(!)" << std::endl;
	}

	// TEST 7 - SHOT MODIFICATION - ROTO-TRANSLATION OF THE SHOT COORDINATES SYSTEM
	if (test8(shot1, shot2, p1, p2))
	{
		std::cout << "TEST 8 (roto-translation of the Shot coordinates system) - PASSED(!)" << std::endl;
	}
	else
	{
		std::cout << "TEST 8 (roto-translation of the Shot coordinates system) - FAILED(!)" << std::endl;
	}	
	
	// TEST 9 - SHOT MODIFICATION - ROTO-TRANSLATION OF THE SHOT COORDINATES SYSTEM
	if (test9(shot1, shot2, p1, p2))
	{
		std::cout << "TEST 9 (roto-translation of the Shot coordinates system) - PASSED(!)" << std::endl;
	}
	else
	{
		std::cout << "TEST 9 (roto-translation of the Shot coordinates system) - FAILED(!)" << std::endl;
	}

  return 0;
}
