bin_PROGRAMS = 	nxspatcher ply2crd crdview nxsbuilder nxsview \
                nxsedit nxsserver nxstest nxs2ply
lib_LIBRARIES = libnexus.a

LIBS=
SDL_FLAGS = -lSDL -lpthread -lGL -lGLU

libnexus_a_SOURCES =  	../../nexus/file.cpp		\
		     	../../nexus/mfile.cpp		\
			../../nexus/nexus.cpp		\
			../../nexus/patch.cpp 		\
			../../nexus/history.cpp 	\
			../../nexus/borderserver.cpp	\
			../../nexus/normalscone.cpp	

nxspatcher_SOURCES =    ../../nexus/crude.cpp        	\
			../../nexus/file.cpp         	\
			../../nexus/mfile.cpp         	\
			../../nexus/nxspatcher.cpp

ply2crd_SOURCES = 	../../nexus/ply2crude.cpp    	\
			../../nexus/crude.cpp        	\
			../../nexus/file.cpp         	\
			../../nexus/mfile.cpp         	\
		  	../../../wrap/ply/plylib.cpp 	



crdview_SOURCES = 	../../nexus/crudeview.cpp       \
			../../nexus/crude.cpp           \
			../../../wrap/gui/trackball.cpp	\
			../../../wrap/gui/trackmode.cpp \
			../../nexus/mfile.cpp         	\
		 	../../nexus/file.cpp

crdview_LDADD = $(SDL_FLAGS) -lGLEW


nxstest_SOURCES = 	../../nexus/file.cpp		\
			../../nexus/nexus.cpp		\
			../../nexus/patch.cpp 		\
		     	../../nexus/mfile.cpp		\
			../../nexus/history.cpp 	\
			../../nexus/normalscone.cpp	\
			../../nexus/borderserver.cpp	\
			../../nexus/watch.cpp		\
			../../nexus/nxstest.cpp		
nxstest_LDADD = -llzo

nxs2ply_SOURCES =       ../../nexus/nxs2ply.cpp        \
			../../nexus/extraction.cpp 	\
		  	../../../wrap/ply/plylib.cpp

nxs2ply_LDADD = libnexus.a -llzo -lGLEW

nxsbuilder_SOURCES = 	../../nexus/nxsbuilder.cpp 	\
			../../nexus/crude.cpp 		\
		     	../../nexus/file.cpp		\
		     	../../nexus/mfile.cpp		\
			../../nexus/nexus.cpp		\
			../../nexus/patch.cpp 		\
			../../nexus/history.cpp 	\
			../../nexus/borderserver.cpp	\
			../../nexus/normalscone.cpp	\
		     	../../nexus/decimate.cpp	\
		     	../../nexus/nxsalgo.cpp		\
			../../nexus/vchain.cpp		\
			../../nexus/vpartition.cpp 	\
			../../nexus/remapping.cpp 	\
			../../nexus/fragment.cpp	\
			../../nexus/watch.cpp		\
			../../nexus/nxsdispatcher.cpp   \
			../../../wrap/ply/plylib.cpp 	

nxsbuilder_LDADD = -llzo -lptypes -lGLEW /usr/local/lib/libANN.a 

nxsserver_SOURCES =     ../../nexus/nxsserver.cpp	\
			../../nexus/fragment.cpp	\
			../../nexus/vpartition.cpp 	\
		     	../../nexus/normalscone.cpp	\
		     	../../nexus/decimate.cpp	

nxsserver_LDADD = -lptypes /usr/local/lib/libANN.a 


nxsview_SOURCES = 	../../nexus/nexusview.cpp 	\
			../../nexus/nexus.cpp 		\
			../../nexus/history.cpp		\
			../../nexus/extraction.cpp 	\
			../../nexus/preload.cpp 	\
			../../nexus/borderserver.cpp	\
			../../nexus/normalscone.cpp	\
			../../nexus/file.cpp 		\
		     	../../nexus/mfile.cpp		\
			../../nexus/nexusmt.cpp 	\
			../../nexus/patch.cpp 		\
		  	../../nexus/watch.cpp		\
		  	../../../wrap/gui/trackball.cpp	\
			../../../wrap/gui/trackmode.cpp 	

nxsview_LDADD = $(SDL_FLAGS) -lGLEW -lglut -llzo -lptypes



nxsedit_SOURCES = 	../../nexus/nxsedit.cpp 	\
			../../nexus/watch.cpp		\
			../../nexus/nexus.cpp 		\
			../../nexus/history.cpp		\
			../../nexus/patch.cpp 		\
			../../nexus/normalscone.cpp	\
			../../nexus/nxsalgo.cpp 	\
			../../nexus/vpartition.cpp      \
			../../nexus/strip.cpp	\
			../../../wrap/strip/tristrip.cpp\
			../../nexus/tristripper/tri_stripper.cpp \
		 	../../nexus/file.cpp 		\
		     	../../nexus/mfile.cpp		\
		  	../../../wrap/ply/plylib.cpp 	\
			../../nexus/borderserver.cpp

nxsedit_LDADD = -llzo -lGLEW /usr/local/lib/libANN.a 

exploder_SOURCES = 	../../nexus/exploder.cpp 	\
			../../nexus/nexus.cpp 		\
			../../nexus/patch.cpp 		\
		 	../../nexus/file.cpp 		\
		     	../../nexus/mfile.cpp		\
			../../nexus/borderserver.cpp

exploder_LDADD = -llzo -lGLEW


INCLUDES = -I../../../
