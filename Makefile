CXX = g++
CFLAGS = -Wall -g -L.
AR       = ar cr
SRC = src/
BIN = bin/
OBJ = obj/
LIBRARY := ${OBJ}lgzstream.a

all:  mkobj mkbin gzstream.o libgzstream.a sff.o ascii.o util.o Read.o Dictionary.o KMerRoutine.o main.o slopmap 
		
#					
slopmap :   $(OBJ)main.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sff.o $(OBJ)gzstream.o
	$(CXX) $(CFLAGS) -o $(BIN)slopmap $(OBJ)main.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sff.o $(OBJ)gzstream.o -I$(LIBRARY) -Xlinker -lz 
	
main.o :  
	$(CXX) -Wall -g  -c -o $(OBJ)main.o $(SRC)main.cpp 
	
KMerRoutine.o :
	$(CXX) -Wall -g  -c -o $(OBJ)KMerRoutine.o $(SRC)KMerRoutine.cpp
	
Dictionary.o :
	$(CXX) -Wall -g  -c -o $(OBJ)Dictionary.o $(SRC)Dictionary.cpp

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

mkobj :
	rm -rf ${OBJ}
	mkdir ${OBJ}

mkbin :
	rm -rf ${BIN}
	mkdir ${BIN}
	

clean:
	rm -rf ${BIN}
	rm -rf ${OBJ}
