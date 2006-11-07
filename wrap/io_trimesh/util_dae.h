#ifndef __VCGLIB_UTILDAE
#define __VCGLIB_UTILDAE

#include <wrap/io_trimesh/additionalinfo.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/allocate.h>

#include <wrap/io_trimesh/io_mask.h>

#include<QtXml/QDomDocument>
#include<QtCore/QFile>
#include <QtCore/QStringList>

#include<vcg/space/point3.h>
#include<vcg/space/tcoord2.h>
#include<vcg/space/color4.h>

namespace vcg {
namespace tri {
namespace io {
	class InfoDAE : public AdditionalInfo
	{
		public:

		InfoDAE()
		{
			mask	= 0;
			numvert = 0;
			numface = 0;
			doc = NULL;
		}

		~InfoDAE()
		{
			delete doc;
			texturefile.clear();
		}

		QDomDocument* doc;		
		std::vector<QString> texturefile; 
	};

	class AdditionalInfoDAE : public AdditionalInfo
	{
	public: 
		vcg::tri::io::InfoDAE* dae;

		AdditionalInfoDAE()
		:AdditionalInfo()
		{
		}

		~AdditionalInfoDAE()
		{
			delete dae;
		}
	};

	class UtilDAE
	{
	public:
		enum DAEError 
		{
			E_NOERROR,				// 0
			E_CANTOPEN,				// 1
			E_NOGEOMETRYLIBRARY,     // 2 
			E_NOMESH,      // 3
			E_NOVERTEXPOSITION,            // 4
			E_NO3DVERTEXPOSITION,			// 5
			E_NO3DSCENE, // 6
			E_INCOMPATIBLECOLLADA141FORMAT, //7
			E_UNREFERENCEBLEDCOLLADAATTRIBUTE, // 8
			E_NOTRIANGLES
		};

		static const char *ErrorMsg(int error)
		{
			static const char * dae_error_msg[] =
			{
				"No errors",
				"Can't open file",
				"File without a geometry library",
				"There isn't mesh in file",
				"The meshes in file haven't the vertex position attribute",
				"The importer assumes that the OpenMeshType uses a 3D point for the vertex position",
				"There isn't any scene in Collada file",
				"The input file is not compatible with COLLADA 1.41 standard format",
				"Collada file is trying to referece an attribute that is not in the file",
				"This version of Collada Importer support only triangular mesh file"
			};

			if(error>9 || error<0) return "Unknown error";
			else return dae_error_msg[error];
		};
	protected:
		inline static void referenceToANodeAttribute(const QDomNode& n,const QString& attr,QString& url_st)
		{
			url_st = n.toElement().attribute(attr);
			int sz = url_st.size() - 1;
			url_st = url_st.right(sz);
			assert(url_st.size() != 0);
		}

		inline static QDomNode findNodeBySpecificAttributeValue(const QDomNode& n,const QString& tag,const QString& attrname,const QString& attrvalue)
		{
			QDomNode ndl = n.toElement();
			return findNodeBySpecificAttributeValue((QDomDocument&) ndl,tag,attrname,attrvalue);
		}

		inline static QDomNode findNodeBySpecificAttributeValue(const QDomDocument& n,const QString& tag,const QString& attrname,const QString& attrvalue)
		{
			QDomNodeList ndl = n.elementsByTagName(tag);
			int ndl_size = ndl.size();
			assert(ndl_size != 0);
			int ind = 0;
			while(ind < ndl_size)
			{
				if (ndl.at(ind).toElement().attribute(attrname) == attrvalue)
					return ndl.at(ind);
				++ind;
			}
			return QDomNode();
		}

		inline static bool isThereTag(const QDomNode& n,const QString& tagname)
		{
			QDomNode ndl = n.toElement();
			return isThereTag((QDomDocument&) n,tagname);
		}

		inline static bool isThereTag(const QDomDocument& n,const QString& tagname)
		{
			return ((n.toElement().elementsByTagName(tagname).size() > 0)? true : false);
		}


		inline static QDomNode attributeSourcePerSimplex(const QDomNode& n,const QDomDocument& startpoint,const QString& sem)
		{
			QDomNodeList vertattr = n.toElement().elementsByTagName("input");
			for(int ind = 0;ind < vertattr.size();++ind)
			{
				if (vertattr.at(ind).toElement().attribute("semantic") == sem)
				{
					QString url; 
					referenceToANodeAttribute(vertattr.at(ind),"source",url);
					return findNodeBySpecificAttributeValue(startpoint,"source","id",url);
				}
			}
			return QDomNode();
		}

		inline static void valueStringList(QStringList& res,const QDomNode& srcnode,const QString& tag) 
		{
			QDomNodeList list = srcnode.toElement().elementsByTagName(tag);
			int list_size = list.size();
			assert(list_size == 1);
			QString nd = list.at(0).firstChild().nodeValue();
			res = nd.split(" ");
			if (res.last() == "")
				res.removeLast();
		
		}

		inline static bool removeChildNodeList(QDomNodeList& nodelst,const QString& tag = "", const QString& attribname = "", const QString& attribvalue = "")
		{
			for(int jj = 0;jj < nodelst.size();++jj)
			{
				removeChildNode(nodelst.at(jj),tag,attribname,attribvalue);
			}
			return true;
		}

	/*	inline static bool removeChildNode(QDomNode& node,const QString& tag = "", const QString& attribname = "", const QString& attribvalue = "")
		{
			return removeChildNode((QDomDocument&) node.toElement(),tag,attribname,attribvalue);
		}*/

		//inline static bool removeChildNode(QDomDocument& node,const QString& tag = "", const QString& attribname = "", const QString& attribvalue = "")
		inline static bool removeChildNode(QDomNode node,const QString& tag = "", const QString& attribname = "", const QString& attribvalue = "")
		{
			QDomNodeList clst = node.childNodes();
			for(int ii = 0;ii < clst.size();++ii)
			{
				QDomNode oldchild = node.childNodes().at(ii); 
				if (tag != "")
				{
					if ((attribname != "") && (attribvalue != ""))
					{
						if (clst.at(ii).toElement().attribute(attribname) == attribvalue)
							node.removeChild(oldchild);
					}
					else if (clst.at(ii).nodeName() == tag) 
							node.removeChild(oldchild);
				}
				else node.removeChild(oldchild);
			}
			return true;
		}

		static void ParseRotationMatrix(vcg::Matrix44f& m,const std::vector<QDomNode>& t)
		{
			vcg::Matrix44f tmp;
			tmp.SetIdentity();
			for(unsigned int ii = 0;ii < t.size();++ii)
			{
				QString rt = t[ii].firstChild().nodeValue();
				QStringList rtl = rt.split(" ");
				if (rtl.last() == "") rtl.removeLast();
				assert(rtl.size() == 4);
				tmp.SetRotate(rtl.at(3).toFloat(),vcg::Point3f(rtl.at(0).toFloat(),rtl.at(1).toFloat(),rtl.at(2).toFloat()));
				tmp *= tmp;
			}
			m = m * tmp;
		}

		static void AddTranslation(vcg::Matrix44f& m,const QDomNode& t)
		{
			QDomNode tr = t.firstChild();
			QString coord = tr.nodeValue();
			QStringList coordlist = coord.split(" ");
			if (coordlist.last() == "") 
				coordlist.removeLast();
			assert(coordlist.size() == 3);
			m[0][0] = 1.0f;
			m[1][1] = 1.0f;
			m[2][2] = 1.0f;
			m[3][3] = 1.0f;
			m[0][3] = coordlist.at(0).toFloat();
			m[1][3] = coordlist.at(1).toFloat();
			m[2][3] = coordlist.at(2).toFloat();
		}

		static void TransfMatrix(const QDomNode& parentnode,const QDomNode& presentnode,vcg::Matrix44f& m)
		{
			if (presentnode == parentnode) return;
			else
			{
				QDomNode par = presentnode.parentNode();
				std::vector<QDomNode> rotlist;
				QDomNode trans;
				for(int ch = 0;ch < par.childNodes().size();++ch)
				{
					if (par.childNodes().at(ch).nodeName() == "rotate")
						rotlist.push_back(par.childNodes().at(ch));
					else if (par.childNodes().at(ch).nodeName() == "translate")
						 {
							trans = par.childNodes().at(ch);
					     }		
				}
				vcg::Matrix44f tmp;
				tmp.SetIdentity();
				if (!trans.isNull()) AddTranslation(tmp,trans);
				ParseRotationMatrix(tmp,rotlist);
				m = m * tmp;
				TransfMatrix(parentnode,par,m);
			}
		}
	};
}
}
}

#endif