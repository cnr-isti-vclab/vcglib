// standard headers
#include <stdio.h>

// stl headers
#include <vector>

// vcg headers
#include<vcg/simplex/vertex/vertex.h>
#include<vcg/simplex/face/with/rtfmfn.h>
#include<vcg/simplex/face/distance.h>
#include<vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/create/platonic.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/space/intersection3.h>

#include <vcg/space/index/aabb_binary_tree/aabb_binary_tree.h>

typedef float AScalarType;

class AEdge;
class AFace;
class AVertex   : public vcg::Vertex< AScalarType, AEdge, AFace > { };
class AFace     : public vcg::FaceRTFMFN< AVertex, AEdge, AFace > { };
class AMesh     : public vcg::tri::TriMesh< std::vector<AVertex>, std::vector<AFace> > { };

typedef vcg::AABBBinaryTreeIndex<AFace, AScalarType, vcg::EmptyClass> AIndex;

static AMesh gMesh;
static AIndex gIndex;

static void CreateMesh(void) {
	vcg::tri::Dodecahedron<AMesh>(gMesh);

	vcg::tri::UpdateFlags<AMesh>::Clear(gMesh);
	vcg::tri::UpdateNormals<AMesh>::PerVertexNormalized(gMesh);
	vcg::tri::UpdateEdges<AMesh>::Set(gMesh);
}

static void SetIndex(void) {
	gIndex.Set(gMesh.face.begin(), gMesh.face.end());
}

static void TestClosest(void) {
	vcg::face::PointDistanceFunctor getPtDist;
	const AIndex::ScalarType x = 0;
	const AIndex::CoordType queryPoint((AIndex::ScalarType)0, (AIndex::ScalarType)0, (AIndex::ScalarType)0);
	const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();

	AIndex::ObjPtr closestFace;
	AIndex::ScalarType closestDist;
	AIndex::CoordType closestPoint;

	closestFace = gIndex.GetClosest(getPtDist, vcg::EmptyClass(), queryPoint, maxDist, closestDist, closestPoint);

	printf("GetClosest Test:\n");

	if (closestFace != 0) {
		printf("\tface     : 0x%p\n", closestFace);
		printf("\tdistance : %f\n", closestDist);
		printf("\tpoint    : [%f, %f, %f]\n", closestPoint[0], closestPoint[1], closestPoint[2]);
	}
	else {
		printf("\tno object found (index is probably empty).\n");
	}
}

static void TestKClosest(void) {
	vcg::face::PointDistanceFunctor getPtDist;
	const unsigned int k = 10;
	const AIndex::CoordType queryPoint((AIndex::ScalarType)0, (AIndex::ScalarType)0, (AIndex::ScalarType)0);
	const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();

	std::vector<AIndex::ObjPtr> closestObjects;
	std::vector<AIndex::ScalarType> closestDistances;
	std::vector<AIndex::CoordType> closestPoints;

	unsigned int rk = gIndex.GetKClosest(getPtDist, vcg::EmptyClass(), k, queryPoint, maxDist, closestObjects, closestDistances, closestPoints);

	printf("GetKClosest Test:\n");
	printf("\tfound %d objects\n", rk);
}

static void TestRay(void) {
	const bool TEST_BACK_FACES = true;

	vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
	const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
	const AIndex::CoordType rayOrigin((AIndex::ScalarType)0, (AIndex::ScalarType)0, (AIndex::ScalarType)0);
	const AIndex::CoordType rayDirection((AIndex::ScalarType)1, (AIndex::ScalarType)0, (AIndex::ScalarType)0);
	const vcg::Ray3<AIndex::ScalarType, false> ray(rayOrigin, rayDirection);

	AIndex::ObjPtr isectFace;
	AIndex::ScalarType rayT;
	AIndex::CoordType isectPt;

	isectFace = gIndex.DoRay(rayIntersector, vcg::EmptyClass(), ray, maxDist, rayT);

	printf("DoRay Test:\n");
	if (isectFace != 0) {
		printf("\tface  : 0x%p\n", isectFace);
		printf("\tray t : %f\n", rayT);
	}
	else {
		printf("\tno object found (index is probably empty).\n");
	}
}

int main (int argc, char ** argv) {
	CreateMesh();

	SetIndex();

	printf("Spatial Index Tests\n");
	printf("---\n");
	TestClosest();
	printf("---\n");
	TestKClosest();
	printf("---\n");
	TestRay();
	printf("---\n");

	return (0);
}
