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

#ifndef VCG_TANGENT_FIELD_OPERATORS
#define VCG_TANGENT_FIELD_OPERATORS

namespace vcg {
	namespace tri{
		
		template <class MeshType>
		class CrossField
		{
			typedef typename MeshType::FaceType FaceType;
			typedef typename MeshType::VertexType VertexType;
			typedef typename MeshType::CoordType CoordType;
			typedef typename MeshType::ScalarType ScalarType;
			typedef typename MeshType::PerFaceAttributeHandle<CoordType> PerFaceAttributeHandle;

		private:
			static ScalarType Sign(ScalarType a){return (ScalarType)((a>0)?+1:-1);}

		public:

			///fird a tranformation matrix to transform 
			///the 3D space to 2D tangent space specified 
			///by the cross field (where Z=0)
			static vcg::Matrix33<ScalarType> TransformationMatrix(FaceType &f)
			{
				typedef typename FaceType::CoordType CoordType;
				typedef typename FaceType::ScalarType ScalarType;

				bool CrossDir0 = vcg::tri::HasPerFaceAttribute(mesh,"CrossDir0");
				assert(CrossDir0);
				Fh0= vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(mesh,std::string("CrossDir0"));

				///transform to 3d
				CoordType axis0=Fh0[&f];
				CoordType axis1=axis0^axis2;
				CoordType axis2=f.N();

				vcg::Matrix33<ScalarType> Trans;

				///it must have right orientation cause of normal
				Trans[0][0]=axis0[0];
				Trans[0][1]=axis0[1];
				Trans[0][2]=axis0[2];
				Trans[1][0]=axis1[0];
				Trans[1][1]=axis1[1];
				Trans[1][2]=axis1[2];
				Trans[2][0]=axis2[0];
				Trans[2][1]=axis2[1];
				Trans[2][2]=axis2[2];

				/////then find the inverse 
				//f.InvTrans=Inverse(f.Trans);
			}

			///transform a given angle from UV (wrt the cross field)
			///to a 3D direction
			static CoordType AngleToVect(const FaceType &f,const ScalarType &angle)
			{
				///find 2D vector
				vcg::Point2<ScalarType> axis2D=vcg::Point2<ScalarType>(cos(angle),sin(angle));
				CoordType axis3D=CoordType(axis2D.X(),axis2D.Y(),0);
				vcg::Matrix33<ScalarType> Trans=TransformationMatrix(f);
				vcg::Matrix33<ScalarType> InvTrans=Inverse(Trans);
				///then transform
				return (InvTrans*axis3D);
			}

			///find an angle with respect to a given face by a given vector
			///in 3D space, it must be projected and normalized with respect to face's normal
			static ScalarType VectToAngle(const FaceType &f,
				const CoordType &vect3D)
			{
				vcg::Matrix33<ScalarType> Trans=TransformationMatrix(f);
				///trensform the vector to the reference frame by rotating it
				CoordType vect_transf=Trans*vect3D;

				///then put to zero to the Z coordinate
				vcg::Point2<ScalarType> axis2D=vcg::Point2<ScalarType>(vect_transf.X(),vect_transf.Y());
				axis2D.Normalize();

				///then find the angle with respact to axis 0
				ScalarType alpha=atan2(axis2D.Y(),axis2D.X());	////to sum up M_PI?
				if (alpha<0)
					alpha=(2*M_PI+alpha);
				if (alpha<0)
					alpha=0;
				return alpha;
			}

			///return the direction of the cross field in 3D
			static void CrossVector(MeshType &mesh,
				const FaceType &f,
				CoordType axis[4])
			{
				bool CrossDir0 = vcg::tri::HasPerFaceAttribute(mesh,"CrossDir0");
				assert(CrossDir0);
				MeshType::PerFaceAttributeHandle<CoordType> Fh0= 
					vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(mesh,std::string("CrossDir0"));
				axis[0]=Fh0[&f];
				axis[1]=f.cN()^axis[0];
				axis[2]=-axis[0];
				axis[3]=-axis[1];
			}

			///return a specific direction given an integer 0..3
			///considering the reference direction of the cross field
			static CoordType CrossVector(MeshType &mesh,
										 const FaceType &f,
										 const int &index)
			{
				assert((index>=0)&&(index<4));
				CoordType axis[4];
				CrossVector(mesh,f,axis);
				return axis[index];
			}

			///rotate a given vector from a face to another
			///vector is expressend in 3d coordinates
			CoordType Rotate(const FaceType &f0,const FaceType &f1,const CoordType &dir3D)
			{
				CoordType N0=f0.cN();
				CoordType N1=f1.cN();

				///find the rotation matrix that maps between normals
				vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(N0,N1);
				CoordType rotated=rotation*dir3D;
				return rotated;
			}

			// returns the 90 deg rotation of a (around n) most similar to target b
			static CoordType K_PI(const CoordType &a, const CoordType &b, const CoordType &n)
			{
				CoordType c = (a^n).normalized();
				ScalarType scorea = a*b;
				ScalarType scorec = c*b;
				if (fabs(scorea)>=fabs(scorec)) return a*Sign(scorea); else return c*Sign(scorec);
			}

			///interpolate cross field with barycentric coordinates
			static CoordType InterpolateCrossField(const CoordType &t0,
													const CoordType &t1,
													const CoordType &t2,
													const CoordType &n,
													const CoordType &bary)
			{
				CoordType trans0=t0;
				CoordType trans1=K_PI(t1,t0,n);
				CoordType trans2=K_PI(t2,t0,n);
				CoordType sum = trans0*bary.X() + trans1 * bary.Y() + trans2 * bary.Z();
				return sum;
			}

			///interpolate cross field with barycentric coordinates using normalized weights
			static typename typename CoordType InterpolateCrossField(const std::vector<CoordType> &TangVect,
				const std::vector<ScalarType> &Weight,
				const std::vector<CoordType> &Norms,
				const typename CoordType &BaseNorm,
				const typename CoordType &BaseDir)
			{
				typedef typename FaceType::CoordType CoordType;
				typedef typename FaceType::ScalarType ScalarType;

				CoordType sum = CoordType(0,0,0);
				for (int i=0;i<TangVect.size();i++)
				{
					CoordType N1=Norms[i];
					///find the rotation matrix that maps between normals
					vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(N1,BaseNorm);
					CoordType rotated=rotation*TangVect[i];
					CoordType Tdir=K_PI(rotated,BaseDir,BaseNorm);
					Tdir.Normalize();
					sum+=(Tdir*Weight[i]);
				}
				sum.Normalize();
				return sum;
			}

			/*///interpolate cross field with barycentric coordinates
			template <class FaceType>
			typename FaceType::CoordType InterpolateCrossField(const typename FaceType::CoordType &t0,
				const typename FaceType::CoordType &t1,
				const typename FaceType::CoordType &n,
				const typename FaceType::ScalarType &weight)
			{
				CoordType trans0=t0;
				CoordType trans1=K_PI(t1,t0,n);
				CoordType sum = t0*weight + MyCross::V( t1, t0, n ) * (1.0-weight);
				return sum;
			}*/
			
			
			///return the difference of two cross field, values between [0,0.5]
			template <class FaceType>
			typename FaceType::ScalarType DifferenceCrossField(const typename FaceType::CoordType &t0,
				const typename FaceType::CoordType &t1,
				const typename FaceType::CoordType &n)
			{
				CoordType trans0=t0;
				CoordType trans1=K_PI(t1,t0,n);
				ScalarType diff = 1-fabs(trans0*trans1);
				return diff;
			}

			///compute the mismatch between 2 faces
			int MissMatch(const FaceType &f0,const FaceType &f1)
			{
				CoordType dir0=CrossVector(f0,0);
				CoordType dir1=CrossVector(f1,0);

				CoordType dir1Rot=Rotate(f1,f0,dir1);
				dir1Rot.Normalize();

				ScalarType angle_diff=VectToAngle(f0,dir1Rot);

				ScalarType step=M_PI/2.0;
				int i=(int)floor((angle_diff/step)+0.5);
				if (i>=0)
					k=i%4;
				else
					k=(-(3*i))%4;
				return k;
			}
			
			///this function return true if a 
			///given vertex is a singular vertex by 
			///moving around i n a roder wai and accounting for 
			///missmatches.. it requires VF topology
			template <class VertexType>
			bool IsSingular(VertexType &v)
			{
				typedef typename VertexType::FaceType FaceType;
				///check that is on border..
				if (v.IsB())
					return false;

				///get first face sharing the edge
				FaceType *f_init=v.VFp();
				int edge_init=v.VFi(); 

				int missmatch=0;
				///and initialize the pos
				vcg::face::Pos<FaceType> VFI(f_init,edge_init);
				bool complete_turn=false;
				do  
				{
					FaceType *curr_f=VFI.F();
					int curr_edge=VFI.E();

					///assert that is not a border edge
					assert(curr_f->FFp(curr_edge)!=curr_f);

					///find the current missmatch
					FaceType *next_f=curr_f->FFp(curr_edge);
					missmatch+=MissMatch(next_f);

					///continue moving 
					VFI.FlipF();
					VFI.FlipE();

					FaceType *next_f=VFI.F();

					///test if I've finiseh with the face exploration
					complete_turn=(next_f==f_init);
					/// or if I've just crossed a mismatch
				}while (!complete_turn);
				return((missmatch%4)!=0);
			}
			
			/*void GetWeights(TriMeshType *mesh,
					const std::vector<FaceType*> &faces,
					std::vector<ScalarType> &weights)
			{
				weights.clear();
				MeshType::PerFaceAttributeHandle<ScalarType> Fh0 FHRVal=
				vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<ScalarType>(mesh,std::string("CrossVal"));

				for (int i=0;i<faces.size();i++)
					weights.push_back(FHRVal[faces[i]]);
			}

			void GetTangDir(TriMeshType *mesh,
						 const std::vector<FaceType*> &faces,
						 std::vector<CoordType> &dir0,
						 std::vector<CoordType> &dir1)
			{
				dir0.clear();
				dir1.clear();
				MeshType::PerFaceAttributeHandle<CoordType> FHDir0=vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(test_mesh,std::string("CrossDir0"));
				MeshType::PerFaceAttributeHandle<CoordType> FHDir1=vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(test_mesh,std::string("CrossDir1"));
				for (int i=0;i<faces.size();i++)
				{
					dir0.push_back(FHDir0[faces[i]]);
					dir1.push_back(FHDir1[faces[i]]);
				}
			}*/

		};///end class
	} //End Namespace Tri
} // End Namespace vcg
#endif