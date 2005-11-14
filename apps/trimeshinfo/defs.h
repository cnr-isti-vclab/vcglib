#ifndef _DEFS_H
#define _DEFS_H
// -----------------------------------------------------------------------------------------------

// error messages
#define MSG_ERR_N_ARGS                  \
  "\nUsage:  "\
  "trimeshinfo filename [opt]\n"\
  "where opt can be:\n"\
  "     -x[y|n]  Enable or not XML info dumping (default no)\n"\
	"     -a[y|n]  Enable or not ascii info dumping (default yes)\n"\
  "     -o<filename.ext>  Save the input mesh in the specified format\n"
			
#define MSG_ERR_MESH_LOAD               "error loading the input meshes.\n"
#define MSG_ERR_INVALID_OPTION          "unable to parse option '%s'\n"

#endif
// -----------------------------------------------------------------------------------------------
