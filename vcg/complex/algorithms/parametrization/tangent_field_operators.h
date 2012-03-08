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

#include <vcg/math/histogram.h>
#include <vcg/complex/algorithms/update/curvature.h>

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

			static void NormalizePerVertImportanceVal(MeshType &mesh)
			{
				vcg::Distribution<ScalarType> Distr;
				for (int i=0;i<mesh.vert.size();i++)
				{
					VertexType *v=&mesh.vert[i];
					if (v->IsD())continue;
					Distr.Add(fabs(v->K2()));
					Distr.Add(fabs(v->K1()));
					/*if (fabs(v->K1())>MaxK)
						MaxK=fabs(fabs(v->K1()));
					if (fabs(v->K2())>MaxK)
						MaxK=fabs(fabs(v->K2()));*/
				}
				ScalarType perc=Distr.Percentile(.99);
				for (int i=0;i<mesh.vert.size();i++)
				{
					VertexType *v=&mesh.vert[i];
					if (v->IsD())continue;
					ScalarType val;
					val=fabs(v->K1());
					if (val>perc)
						val=perc;
					else
						val/=perc;

					v->K1()=val;
					val=(v->K2());
					if (val>perc)
						val=perc;
					else
						val/=perc;
					v->K2()=val;
				}
			}
		
			//static void NormalizePerVertImportanceVal(MeshType &mesh)
			//{
			//	//vcg::Distribution<ScalarType> Distr;
			//	ScalarType MaxK=0;
			//	ScalarType MinK=0;
			//	for (int i=0;i<mesh.vert.size();i++)
			//	{
			//		VertexType *v=&mesh.vert[i];
			//		if (v->IsD())continue;
			//		Distr.Add((v->K2()));
			//		Distr.Add((v->K1()));
			//		if (fabs(v->K1())>MaxK)
			//			MaxK=fabs(fabs(v->K1()))
			//		if (fabs(v->K2())>MaxK)
			//			MaxK=fabs(fabs(v->K2()));*/
			//	}
			//	ScalarType perc0=Distr.Percentile(.1);
			//	ScalarType perc1=Distr.Percentile(.9);
			//	ScalarType val=perc0-perc1
			//	for (int i=0;i<mesh.vert.size();i++)
			//	{
			//		VertexType *v=&mesh.vert[i];
			//		if (v->IsD())continue;
			//		ScalarType val;
			//		val=(v->K1());
			//		if (val<perc0)
			//			val=perc0;
			//		else
			//			val/=perc;
			//		v->K1()=val;
			//		val=(v->K2());
			//		if (val>perc)
			//			val=perc;
			//		else
			//			val/=perc;
			//		v->K2()=val;
			//	}
			//}

		public:
			static void AddCrossAttributesIfNeeded(MeshType &mesh,
				PerFaceAttributeHandle &_FHDir0,
				PerFaceAttributeHandle &_FHDir1)
				//PerFaceScalarAttributeHandle &_FHVal)
			{
				bool CrossDir0 = vcg::tri::HasPerFaceAttribute(mesh,"CrossDir0");
				bool CrossDir1 = vcg::tri::HasPerFaceAttribute(mesh,"CrossDir1");
				bool CrossVal = vcg::tri::HasPerFaceAttribute(mesh,"CrossVal");

				if (!CrossDir0)
					_FHDir0=vcg::tri::Allocator<MeshType>::AddPerFaceAttribute<CoordType>(mesh,std::string("CrossDir0"));
				else
					_FHDir0=vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(mesh,std::string("CrossDir0"));

				if (!CrossDir1)
					_FHDir1=vcg::tri::Allocator<MeshType>::AddPerFaceAttribute<CoordType>(mesh,std::string("CrossDir1"));
				else
					_FHDir1=vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(mesh,std::string("CrossDir1"));


				/*if (!CrossVal)
				_FHVal=vcg::tri::Allocator<MeshType>::AddPerFaceAttribute<ScalarType>(test_mesh,std::string("CrossVal"));
				else
				_FHVal=vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<ScalarType>(test_mesh,std::string("CrossVal"));*/
			}

			static void CopyFaceCrossField(MeshType &left,MeshType &right)
			{

				PerFaceAttributeHandle _FHR0,_FHR1;
				PerFaceScalarAttributeHandle _FHRVal;
				bool CrossDir0 = vcg::tri::HasPerFaceAttribute(right,"CrossDir0");
				bool CrossDir1 = vcg::tri::HasPerFaceAttribute(right,"CrossDir1");
				//bool CrossVal = vcg::tri::HasPerFaceAttribute(right,"CrossVal");
				assert(CrossDir0);
				assert(CrossDir1);
				//assert(CrossVal);

				_FHR0=vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(right,std::string("CrossDir0"));
				_FHR1=vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(right,std::string("CrossDir1"));
				//_FHRVal=vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<ScalarType>(right,std::string("CrossVal"));

				PerFaceAttributeHandle _FHL0,_FHL1;
				//PerFaceScalarAttributeHandle _FHLVal;
				AddCrossAttributesIfNeeded(left,_FHL0,_FHL1);//,_FHLVal);

				assert(left.face.size()==right.face.size());
				for (int i=0;i<left.face.size();i++)
				{
					_FHL0[i]=_FHR0[i];
					_FHL1[i]=_FHR1[i];
					//_FHLVal[i]=_FHRVal[i];
				}
			}

			static void SetVertCrossFromCurvature(MeshType &mesh)
			{
				vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
				vcg::tri::UpdateTopology<MeshType>::VertexFace(mesh);
				vcg::tri::UpdateBounding<MeshType>::Box(mesh);

				//set as selected high curvature value
				vcg::tri::UpdateCurvature<MeshType>::PrincipalDirectionsNormalCycles(mesh);
				NormalizePerVertImportanceVal(mesh);
				///save the curvature value
				std::vector<ScalarType> K1,K2;
				K1.resize(mesh.vert.size());
				K2.resize(mesh.vert.size());
				for (int j=0;j<mesh.vert.size();j++)
				{
					VertexType *v=&mesh.vert[j];
					if(v->IsD())continue;
					K1[j]=v->K1();
					K2[j]=v->K2();
				}
				///then find multiscale curvature directions
				vcg::tri::UpdateCurvature<MeshType>::PrincipalDirectionsPCA(mesh,mesh.bbox.Diag()/200.0);
				///and save back importance val
				for (int j=0;j<mesh.vert.size();j++)
				{
					VertexType *v=&mesh.vert[j];
					if(v->IsD())continue;
					v->K1()=K1[j];
					v->K2()=K2[j];
				}

				///set normal according to curvature
				for (int j=0;j<mesh.vert.size();j++)
				{
					VertexType *v=&mesh.vert[j];
					if(v->IsD())continue;
					CoordType N0=v->N();
					v->N()=v->PD1()^v->PD2();
					v->N().Normalize();
					if (N0*v->N()<0)
						v->N()=-v->N();
				}
			}

			///fird a tranformation matrix to transform 
			///the 3D space to 2D tangent space specified 
			///by the cross field (where Z=0)
			static vcg::Matrix33<ScalarType> TransformationMatrix(MeshType &mesh,const FaceType &f)
			{
				typedef typename FaceType::CoordType CoordType;
				typedef typename FaceType::ScalarType ScalarType;

				bool CrossDir0 = vcg::tri::HasPerFaceAttribute(mesh,"CrossDir0");
				assert(CrossDir0);
				PerFaceAttributeHandle Fh0= vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(mesh,std::string("CrossDir0"));

				///transform to 3d
				CoordType axis0=Fh0[&f];
				CoordType axis1=axis0^f.cN();
				CoordType axis2=f.cN();

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
				return (Trans);
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
			static ScalarType VectToAngle(MeshType &mesh,const FaceType &f,const CoordType &vect3D)
			{
				vcg::Matrix33<ScalarType> Trans=TransformationMatrix(mesh,f);
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
			///given a first direction
			static void CrossVector(const CoordType &dir0,
				const CoordType &norm,
				CoordType axis[4])
			{
				axis[0]=dir0;
				axis[1]=norm^axis[0];
				axis[2]=-axis[0];
				axis[3]=-axis[1];
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
				CoordType dir0=Fh0[&f];
				CrossVector(dir0,f.cN(),axis);
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

			///return the direction of the cross field in 3D
			static void SetCrossVector(MeshType &mesh,
				const FaceType &f,
				CoordType dir0,
				CoordType dir1)
			{
				bool CrossDir0 = vcg::tri::HasPerFaceAttribute(mesh,"CrossDir0");
				assert(CrossDir0);
				bool CrossDir1 = vcg::tri::HasPerFaceAttribute(mesh,"CrossDir1");
				assert(CrossDir1);
				MeshType::PerFaceAttributeHandle<CoordType> Fh0=
					vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(mesh,std::string("CrossDir0"));
				MeshType::PerFaceAttributeHandle<CoordType> Fh1=
					vcg::tri::Allocator<MeshType>::GetPerFaceAttribute<CoordType>(mesh,std::string("CrossDir1"));
				Fh0[f]=dir0;
				Fh1[f]=dir1;
			}

			///rotate a given vector from a face to another
			///vector is expressend in 3d coordinates
			static CoordType Rotate(const FaceType &f0,const FaceType &f1,const CoordType &dir3D)
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
													const CoordType &n0,
													const CoordType &n1,
													const CoordType &n2,
													const CoordType &target_n,
													const CoordType &bary)
			{
				vcg::Matrix33<ScalarType> R0=vcg::RotationMatrix(n0,target_n);
				vcg::Matrix33<ScalarType> R1=vcg::RotationMatrix(n1,target_n);
				vcg::Matrix33<ScalarType> R2=vcg::RotationMatrix(n2,target_n);
				///rotate
				CoordType trans0=R0*t0;
				CoordType trans1=R1*t1;
				CoordType trans2=R2*t2;
				/*CoordType trans0=t0;
				CoordType trans1=t1;
				CoordType trans2=t2;*/
				trans0.Normalize();
				trans1.Normalize();
				trans2.Normalize();
				///k_PI/2 rotation
				trans1=K_PI(trans1,trans0,target_n);
				trans2=K_PI(trans2,trans0,target_n);
				trans1.Normalize();
				trans2.Normalize();

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

			///interpolate cross field with scalar weight
			static typename FaceType::CoordType InterpolateCrossFieldLine(const typename FaceType::CoordType &t0,
																const typename FaceType::CoordType &t1,
																const typename FaceType::CoordType &n0,
																const typename FaceType::CoordType &n1,
																const typename FaceType::CoordType &target_n,
																const typename FaceType::ScalarType &weight)
			{
				vcg::Matrix33<ScalarType> R0=vcg::RotationMatrix(n0,target_n);
				vcg::Matrix33<ScalarType> R1=vcg::RotationMatrix(n1,target_n);
				CoordType trans0=R0*t0;
				CoordType trans1=R1*t1;
				//CoordType trans0=t0;//R0*t0;
				//CoordType trans1=t1;//R1*t1;
				trans0.Normalize();
				trans1.Normalize();
				trans1=K_PI(trans1,trans0,target_n);
				trans1.Normalize();
				CoordType sum = trans0*weight + trans1 * (1.0-weight);
				return sum;
			}


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
			static int MissMatch(MeshType &mesh,
				const FaceType &f0,
				const FaceType &f1)
			{
				CoordType dir0=CrossVector(mesh,f0,0);
				CoordType dir1=CrossVector(mesh,f1,0);

				CoordType dir1Rot=Rotate(f1,f0,dir1);
				dir1Rot.Normalize();

				ScalarType angle_diff=VectToAngle(mesh,f0,dir1Rot);

				ScalarType step=M_PI/2.0;
				int i=(int)floor((angle_diff/step)+0.5);
				int k=0;
				if (i>=0)
					k=i%4;
				else
					k=(-(3*i))%4;
				return k;
			}

			static void SortedFaces(MeshType &mesh,
				VertexType &v,
				std::vector<FaceType*> &faces)
			{
				typedef typename VertexType::FaceType FaceType;

				///check that is not on border..
				assert (!v.IsB());

				///get first face sharing the edge
				FaceType *f_init=v.VFp();
				int edge_init=v.VFi(); 

				///and initialize the pos
				vcg::face::Pos<FaceType> VFI(f_init,edge_init);
				bool complete_turn=false;
				do  
				{
					FaceType *curr_f=VFI.F();
					faces.push_back(curr_f);

					int curr_edge=VFI.E();

					///assert that is not a border edge
					assert(curr_f->FFp(curr_edge)!=curr_f);

					/*///find the current missmatch
					missmatch+=(curr_f,const FaceType &f1);*/

					///continue moving 
					VFI.FlipF();
					VFI.FlipE();

					FaceType *next_f=VFI.F();

					///test if I've finiseh with the face exploration
					complete_turn=(next_f==f_init);
					/// or if I've just crossed a mismatch
				}while (!complete_turn);
			}

			/////this function return true if a 
			/////given vertex is a singular vertex by 
			/////moving around i n a roder wai and accounting for 
			/////missmatches.. it requires VF topology
			/////this function return true if a 
			/////given vertex is a singular vertex by 
			/////moving around i n a roder wai and accounting for 
			/////missmatches
			//static bool IsSingular(MeshType &mesh,
			//						VertexType &v,
			//						int &missmatch)
			//{
			//	typedef typename VertexType::FaceType FaceType;
			//	///check that is on border..
			//	if (v.IsB())
			//		return false;

			//	///get first face sharing the edge
			//	FaceType *f_init=v.VFp();
			//	int edge_init=v.VFi(); 

			//	//int missmatch=0;
			//	missmatch=0;
			//	///and initialize the pos
			//	vcg::face::Pos<FaceType> VFI(f_init,edge_init);
			//	bool complete_turn=false;
			//	do  
			//	{
			//		FaceType *curr_f=VFI.F();
			//		int curr_edge=VFI.E();

			//		///assert that is not a border edge
			//		assert(curr_f->FFp(curr_edge)!=curr_f);

			//		/*///find the current missmatch
			//		missmatch+=(curr_f,const FaceType &f1);*/

			//		///continue moving 
			//		VFI.FlipF();
			//		VFI.FlipE();

			//		FaceType *next_f=VFI.F();
			//		
			//		///find the current missmatch
			//		missmatch+=MissMatch(mesh,*curr_f,*next_f);

			//		///test if I've finiseh with the face exploration
			//		complete_turn=(next_f==f_init);
			//		/// or if I've just crossed a mismatch
			//	}while (!complete_turn);
			//	missmatch=missmatch%4;
			//	return(missmatch!=0);
			//}

			static int SimilarDir(CoordType dir[4],
				CoordType TestD)
			{
				int ret=-1;
				ScalarType maxAcc=-1;
				for (int i=0;i<4;i++)
				{
					ScalarType testAcc=fabs(dir[i]*TestD);
					if (testAcc>maxAcc)
					{
						maxAcc=testAcc;
						ret=i;
					}
				}
				assert(ret!=-1);
				return ret;
			}

			static bool IsSingular(MeshType &mesh,
				VertexType &v,
				int &missmatch)
			{
				typedef typename VertexType::FaceType FaceType;
				///check that is on border..
				if (v.IsB())
					return false;

				std::vector<FaceType*> faces;
				SortedFaces(mesh,v,faces);
				for (int i=0;i<faces.size();i++)
				{
					FaceType *curr_f=faces[i];
					FaceType *next_f=faces[(i+1)%faces.size()];

					///find the current missmatch
					missmatch+=MissMatch(mesh,*curr_f,*next_f);

				}
				missmatch=missmatch%4;
				return(missmatch!=0);
			}

			static bool LoadFIELD(MeshType *mesh,
				const char *path_vfield,
				bool per_vertex=false)
			{
				FILE *f = fopen(path_vfield,"rt");
				if (!f) {
					//if (errorMsg) sprintf(errorMsg,"Cannot Open File :(");
					return false;
				}
				{
					char word[512]; word[0]=0;
					fscanf(f,"%s",word);
					char c=0;
					if (word[0]=='#') {
						// skip comment line
						while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break;
					} else {
						//if (errorMsg) sprintf(errorMsg,"The VField file should start with a comment");
						return false;
					}
					int nnv = -1;
					if (fscanf(f,"%d",&nnv)!=1) {
						// number of vertices not read. Skip another line (ffield file?) and try again.
						while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break; // skip
						fscanf(f,"%d",&nnv);
					}
					int targetnum=mesh->fn;
					if (per_vertex)
						targetnum=mesh->vn;
					if (nnv != (int)targetnum) 
					{
						//if (errorMsg) sprintf(errorMsg,"Wrong element number. Found: %d. Expected: %d.",nnv,mesh->vn);
						return false;
					}
					while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break; // skip
					// skip strange string line
					while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break;
					for (int i=0; i<nnv; i++){
						vcg::Point3d u,v;
						int a,b;
						if (fscanf(f,
							"%d %d %lf %lf %lf %lf %lf %lf",
							&a,&b,
							&(v.X()),&(v.Y()),&(v.Z()),
							&(u.X()),&(u.Y()),&(u.Z())
							)!=8) {
								//if (errorMsg) sprintf(errorMsg,"Format error reading vertex n. %d",i);
								return false;
						}
						//node[i]->TF().Import(u);
						if (per_vertex)
						{
							mesh->vert[i].PD1()=u;
							mesh->vert[i].PD2()=v;
						}
						else
						{
							FaceType *f=&mesh->face[i];
							SetCrossVector(*mesh,*f,u,v);
						}
					}
				}
				fclose(f);
				return true;
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