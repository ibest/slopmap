CXX = g++
CFLAGS = -Wall -g -L.
AR       = ar cr
SRC = src/
SRC_STL = src_stl/
BIN = bin/
OBJ = obj/
LIBRARY := ${OBJ}lgzstream.a

all:  mkobj mkbin gzstream.o libgzstream.a sff.o ascii.o util.o Read.o dnautil.o Dictionary.o KMerRoutine.o main.o slopmap 
		
#					
slopmap :   $(OBJ)main.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)dnautil.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sff.o $(OBJ)gzstream.o
	$(CXX) $(CFLAGS) -o $(BIN)slopmap $(OBJ)main.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)dnautil.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sff.o $(OBJ)gzstream.o -I$(LIBRARY) -Xlinker -lz
	
main.o :  
	$(CXX) -Wall -g  -c -o $(OBJ)main.o $(SRC)main.cpp 
	
KMerRoutine.o :
	$(CXX) -Wall -g  -c -o $(OBJ)KMerRoutine.o $(SRC)KMerRoutine.cpp
	
Dictionary.o :
	$(CXX) -Wall -g  -c -o $(OBJ)Dictionary.o $(SRC)Dictionary.cpp
	
dnautil.o :  
	$(CXX) -Wall -g  -c -o $(OBJ)dnautil.o $(SRC)dnautil.cpp 

Read.o :
	$(CXX) -Wall -g  -c -o $(OBJ)Read.o $(SRC)Read.cpp
	
util.o :
	$(CXX) -Wall -g  -c -o $(OBJ)util.o $(SRC)util.cpp
	
ascii.o :
	$(CXX) -Wall -g  -c -o $(OBJ)ascii.o $(SRC)ascii.cpp
	
sff.o: $(SRC)sff.h $(SRC)sff.c
	g++ -g -I $(SRC) -c -o $(OBJ)sff.o $(SRC)sff.c
	
	
libgzstream.a : $(OBJ)gzstream.o $(SRC)gzstream.h
	${AR} $(OBJ)libgzstream.a $(OBJ)gzstream.o
	
gzstream.o : $(SRC)gzstream.C $(SRC)gzstream.h
	$(CXX) -I $(SRC) -O -Wall -c -o $(OBJ)gzstream.o $(SRC)gzstream.C 


#Target "stl" - only standard template libraty is used. Needed to build the SlopMap on some machines
 
stl : mkobj mkbin gzstream_stl libgzstream_stl sff_stl ascii_stl util_stl Read_stl dnautil_stl Dictionary_stl KMerRoutine_stl main_stl slopmap_stl 
		
#					
slopmap_stl :   $(OBJ)main.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)dnautil.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sff.o $(OBJ)gzstream.o
	$(CXX) $(CFLAGS) -o $(BIN)slopmap $(OBJ)main.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)dnautil.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sff.o $(OBJ)gzstream.o -I$(LIBRARY) -Xlinker -lz
	
main_stl :  
	$(CXX) -Wall -g  -c -o $(OBJ)main.o $(SRC_STL)main.cpp 
	
KMerRoutine_stl :
	$(CXX) -Wall -g  -c -o $(OBJ)KMerRoutine.o $(SRC_STL)KMerRoutine.cpp
	
Dictionary_stl :
	$(CXX) -Wall -g  -c -o $(OBJ)Dictionary.o $(SRC_STL)Dictionary.cpp
	
dnautil_stl :  
	$(CXX) -Wall -g  -c -o $(OBJ)dnautil.o $(SRC_STL)dnautil.cpp 

Read_stl :
	$(CXX) -Wall -g  -c -o $(OBJ)Read.o $(SRC_STL)Read.cpp
	
util_stl :
	$(CXX) -Wall -g  -c -o $(OBJ)util.o $(SRC_STL)util.cpp
	
ascii_stl :
	$(CXX) -Wall -g  -c -o $(OBJ)ascii.o $(SRC_STL)ascii.cpp
	
sff_stl : $(SRC_STL)sff.h $(SRC_STL)sff.c
	g++ -g -I $(SRC_STL) -c -o $(OBJ)sff.o $(SRC_STL)sff.c
	
	
libgzstream_stl : $(OBJ)gzstream.o $(SRC_STL)gzstream.h
	${AR} $(OBJ)libgzstream.a $(OBJ)gzstream.o
	
gzstream_stl : $(SRC_STL)gzstream.C $(SRC_STL)gzstream.h
	$(CXX) -I $(SRC_STL) -O -Wall -c -o $(OBJ)gzstream.o $(SRC_STL)gzstream.C 

#-------------------------------------------------------------------------------

mkobj :
	rm -rf ${OBJ}
	mkdir ${OBJ}

mkbin :
	rm -rf ${BIN}
	mkdir ${BIN}
	

	

clean:
	rm -rf ${BIN}
	rm -rf ${OBJ}
