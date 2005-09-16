CC	 	= g++
CFLAGS = -Wall -pedantic
DEFINED = __GNUC
ROOT=../../../..


app: plylib.o $(ROOT)/apps/metro/metro.cpp
	$(CC) $(CFLAGS) -D$(DEFINED) -I$(ROOT) -o metro $(ROOT)/apps/metro/metro.cpp plylib.o
	
plylib.o: $(ROOT)/wrap/ply/plylib.cpp
	$(CC) -c -D$(DEFINED) -I$(ROOT) $(ROOT)/wrap/ply/plylib.cpp

	
