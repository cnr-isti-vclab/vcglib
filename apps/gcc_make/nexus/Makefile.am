bin_PROGRAMS = 	nxspatcher ply2crd crdview voronoinxs nxsview nxsedit 

LIBS=
SDL_FLAGS = -lSDL -lpthread -lGL -lGLU

nxspatcher_SOURCES =    ../../nexus/crude.cpp        	\
			../../nexus/file.cpp         	\
			../../nexus/nxspatcher.cpp

ply2crd_SOURCES = 	../../nexus/ply2crude.cpp    	\
			../../nexus/crude.cpp        	\
			../../nexus/file.cpp         	\
		  	../../../wrap/ply/plylib.cpp 	



crdview_SOURCES = 	../../nexus/crudeview.cpp       \
			../../nexus/crude.cpp           \
			../../nexus/vert_remap.cpp      \
			../../nexus/mfhash.cpp 		\
			../../../wrap/gui/trackball.cpp	\
			../../../wrap/gui/trackmode.cpp \
		 	../../nexus/file.cpp

crdview_LDADD = $(SDL_FLAGS) -lGLEW


#test_SOURCES = 		../../nexus/file.cpp		\
#			../../nexus/nexus.cpp		\
#			../../nexus/patch.cpp 		\
#		    	../../nexus/patchserver.cpp	\
#			../../nexus/borderserver.cpp	\
#			../../nexus/test.cpp		
#test_LDADD = -llzo

voronoinxs_SOURCES = 	../../nexus/voronoinxs.cpp 	\
			../../nexus/crude.cpp 		\
		     	../../nexus/file.cpp		\
		     	../../nexus/mfile.cpp		\
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
			../../nexus/fragment.cpp	\
			../../nexus/watch.cpp		\
			../../../wrap/ply/plylib.cpp 	

voronoinxs_LDADD = -llzo -lptypes -lGLEW


nxsview_SOURCES = 	../../nexus/nexusview.cpp 	\
			../../nexus/nexus.cpp 		\
		    	../../nexus/patchserver.cpp	\
			../../nexus/borderserver.cpp	\
			../../nexus/file.cpp 		\
		     	../../nexus/mfile.cpp		\
			../../nexus/nexusmt.cpp 	\
			../../nexus/patch.cpp 		\
		  	../../../wrap/gui/trackball.cpp	\
		  	../../nexus/watch.cpp	\
			../../../wrap/gui/trackmode.cpp 	

nxsview_LDADD = $(SDL_FLAGS) -lGLEW -lglut -llzo



nxsedit_SOURCES = 	../../nexus/nxsedit.cpp 	\
			../../nexus/watch.cpp		\
			../../nexus/nexus.cpp 		\
			../../nexus/patch.cpp 		\
			../../nexus/nxsalgo.cpp 	\
			../../../wrap/strip/tristrip.cpp\
			../../nexus/tristripper/tri_stripper.cpp \
		 	../../nexus/file.cpp 		\
		     	../../nexus/mfile.cpp		\
		  	../../../wrap/ply/plylib.cpp 	\
			../../nexus/patchserver.cpp 	\
			../../nexus/borderserver.cpp

nxsedit_LDADD = -llzo -lGLEW


INCLUDES = -I../../../
