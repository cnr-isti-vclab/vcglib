bin_PROGRAMS = 	ply2crd crdview voronoinxs nxsview nxsedit


SDL_FLAGS = -lSDL -lpthread -lGL -lGLU


ply2crd_SOURCES = 	../../nexus/ply2crude.cpp    	\
			../../nexus/crude.cpp        	\
			../../nexus/file.cpp         	\
		  	../../../wrap/ply/plylib.cpp 	



crdview_SOURCES = 	../../nexus/crudeview.cpp       \
			../../nexus/crude.cpp           \
			../../nexus/vert_remap.cpp      \
			../../nexus/mfhash.cpp 		\
		 	../../nexus/file.cpp
crdview_LDADD = $(SDL_FLAGS)



voronoinxs_SOURCES = 	../../nexus/voronoinxs.cpp 	\
			../../nexus/crude.cpp 		\
		     	../../nexus/file.cpp		\
			../../nexus/vert_remap.cpp	\
			../../nexus/mfhash.cpp 		\
			../../nexus/nexus.cpp		\
			../../nexus/patch.cpp 		\
		    	../../nexus/patchserver.cpp	\
			../../nexus/borderserver.cpp	\
		     	../../nexus/decimate.cpp	\
			../../nexus/nxsbuild.cpp 	\
			../../nexus/voronoichain.cpp	\
			../../nexus/pvoronoi.cpp 	\
			../../../wrap/ply/plylib.cpp 	



nxsview_SOURCES = 	../../nexus/nexusview.cpp 	\
			../../nexus/nexus.cpp 		\
		    	../../nexus/patchserver.cpp	\
			../../nexus/borderserver.cpp	\
			../../nexus/file.cpp 		\
			../../nexus/nexusmt.cpp 	\
			../../nexus/patch.cpp 		\
		  	../../../wrap/gui/trackball.cpp	\
			../../../wrap/gui/trackmode.cpp 	

nxsview_LDFLAGS = $(SDL_FLAGS) -lGLEW



nxsedit_SOURCES = 	../../nexus/nxsedit.cpp 	\
			../../nexus/nexus.cpp 		\
			../../nexus/patch.cpp 		\
			../../nexus/nxsalgo.cpp 	\
			../../../wrap/strip/tristrip.cpp\
			../../nexus/tristripper/tri_stripper.cpp \
		 	../../nexus/file.cpp 		\
			../../nexus/patchserver.cpp 	\
			../../nexus/borderserver.cpp



INCLUDES = -I../../../
