CC	 	= g++
CFLAGS = -Wall -pedantic
DEFINED = __GNUC
ROOT=../../../..


app: plylib.o $(ROOT)/apps/tri_decimator/main.cpp
	$(CC) $(CFLAGS) -D$(DEFINED) -I$(ROOT) -o tri_decimator $(ROOT)/apps/tri_decimator/main.cpp plylib.o
	
plylib.o: $(ROOT)/wrap/ply/plylib.cpp
	$(CC) -c -D$(DEFINED) -I$(ROOT) $(ROOT)/wrap/ply/plylib.cpp

	
