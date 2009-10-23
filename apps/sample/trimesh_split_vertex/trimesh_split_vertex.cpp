#include <vector>

#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/attribute_seam.h>

#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

/*
  this sample shows how to transfer per wedge attributes from wedges to vertices.
  during the process new vertices could be created.
*/

// source mesh type: per-wedge texture coordinates
class SrcVertex;
class SrcEdge;
class SrcFace;

class SrcVertex : public vcg::VertexSimp2   <SrcVertex, SrcEdge, SrcFace, vcg::vertex::Coord3f> { };
class SrcFace   : public vcg::FaceSimp2     <SrcVertex, SrcEdge, SrcFace, vcg::face::VertexRef, vcg::face::WedgeTexCoord2f> { };
class SrcMesh   : public vcg::tri::TriMesh  <std::vector<SrcVertex>, std::vector<SrcFace> > { };


// destination mesh type: per-vertex texture coordinates
class DstVertex;
class DstEdge;
class DstFace;

class DstVertex : public vcg::VertexSimp2   <DstVertex, DstEdge, DstFace, vcg::vertex::Coord3f, vcg::vertex::TexCoord2f> { };
class DstFace   : public vcg::FaceSimp2     <DstVertex, DstEdge, DstFace, vcg::face::VertexRef> { };
class DstMesh   : public vcg::tri::TriMesh  <std::vector<DstVertex>, std::vector<DstFace> > { };

// extract wedge attributes functor.
// given a source face and a wedge index, this functor extracts all the relevant attributes from the wedge
// and transfer them to the destination vertex.
// source and destination meshes are provided to allow for attribute presence checking (.Is*Enabled()).
inline void ExtractVertex(const SrcMesh & srcMesh, const SrcFace & f, int whichWedge, const DstMesh & dstMesh, DstVertex & v)
{
	(void)srcMesh;
	(void)dstMesh;

	v.P() = f.cP(whichWedge);
	v.T() = f.cWT(whichWedge);
}

// sample compare functor.
// given two destination vertices, this functor tells if they are identical in all relevan attributes.
// source and destination meshes are provided to allow for attribute presence checking (.Is*Enabled()).
inline bool CompareVertex(const DstMesh & m, const DstVertex & vA, const DstVertex & vB)
{
	(void)m;

	return (vA.cT() == vB.cT());
}

// sample copy functor.
// given two destination vertices, this functor is asked to copy all relevan attributes.
// source and destination meshes are provided to allow for attribute presence checking (.Is*Enabled()).
inline void CopyVertex(const DstMesh & m, const DstVertex & vSrc, DstVertex & vDst)
{
	(void)m;

	vDst.P() = vSrc.cP();
	vDst.T() = vSrc.cT();
}

void usage(void)
{
	printf("usage : trimesh_split_vertex <src_ply_file_name> <dst_ply_file_name>\n");
	printf("where : <src_ply_file_name> : source PLY trimesh file name with texture coordinates per wedge\n");
	printf("        <dst_ply_file_name> : destination PLY trimesh file name with texture coordinates per vertex\n");
	printf("exit.\n");
}

int main(int argc, char ** argv)
{
	if (argc != 3)
	{
		usage();
		return -1;
	}

	SrcMesh srcMesh;
	vcg::tri::io::ImporterPLY<SrcMesh>::Open(srcMesh, argv[1]);
	if ((srcMesh.vn <= 0) || (srcMesh.fn <= 0))
	{
		printf("invalid source mesh file.\n");
		return -1;
	}
	printf("source mesh succesfully loaded.\n");

	DstMesh dstMesh;
	vcg::tri::AttributeSeam::SplitVertex(srcMesh, dstMesh, ExtractVertex, CompareVertex, CopyVertex);
	dstMesh.textures = srcMesh.textures;
	if (vcg::tri::io::ExporterPLY<DstMesh>::Save(dstMesh, argv[2], vcg::tri::io::Mask::IOM_VERTCOORD | vcg::tri::io::Mask::IOM_VERTTEXCOORD) != 0)
	{
		printf("cannot save destination mesh file.\n");
		return -1;
	}
	printf("destination mesh succesfully saved.\n");

	printf("\n");
	printf("statistics:\n");
	printf("  input mesh vertices count    : %d\n", srcMesh.vn);
	printf("  input mesh faces count       : %d\n", srcMesh.fn);
	printf("  splitted mesh vertices count : %d\n", dstMesh.vn);
	printf("  splitted mesh faces count    : %d\n", dstMesh.fn);
	printf("\n");

	return 0;
}
