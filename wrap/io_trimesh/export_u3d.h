#ifndef __VCGLIB_EXPORTERU3D
#define __VCGLIB_EXPORTERU3D

#include <string>

namespace vcg {
namespace tri {
namespace io {

template<typename SaveMeshType>
class ExporterU3D
{
private:

#ifdef WIN32
	static void CreateDir()
	{
		system("mkdir __tmp_dir");
	}

	static void RemoveDir()
	{
		system("rmdir /S /Q __tmp_dir");
	}

	static void InvokeConverter(const char* converter_path,const char* idtf_file,const char* u3d_file)
	{
		std::string com(converter_path);
		system(("cd " + com).c_str());
		com = std::string("IDTFConverter.exe -input ") + idtf_file + " -output " + u3d_file;
		system(com.c_str());
	}

#else
	#ifdef LINUX
	static void CreateDir()
	{
		system("echo off mkdir __tmp_dir");
	}

	static void RemoveDir()
	{
		system("echo off rmdir -ignore-fail-on-non-empty __tmp_dir");
	}

	static void InvokeConverter(const char* converter_path,const char* idtf_file,const char* u3d_file)
	{

	}
	#endif
#endif

public:
	
	static void Save(SaveMeshType& m,const char* outputfile,const int mask)
	{
		CreateDir();
		std::string tmp(outputfile);
		tmp = "tmp_dir/" + tmp + ".idtf";
		vcg::tri::io::ExporterIDTF<SaveMeshType>::Save(m,tmp.c_str(),mask);
		InvokeConverter("../../code/lib/U3D/Bin/Win32/Release",tmp.c_str() ,outputfile);
		RemoveDir();
	}
};

}
}
}

#endif