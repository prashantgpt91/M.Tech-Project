SRC = BaseClass.cpp Configuration.cpp INDetector.cpp Garg.cpp Graph.cpp GraphMatcher.cpp Igraph.cpp Matrix.cpp Sampler.cpp OtherSimilarity.cpp

CC = g++
CPPFLAGS = -O2 -W -Wall -g
OBJ = $(SRC:.cpp=.o)
BIN = ../Bins/INDectector_linux
LDLIBS = -lm -lz

.c.o:

all: $(OBJ)
	$(CC) -o $(BIN) $(OBJ) $(LDLIBS) 

clean:
	rm -f $(OBJ)

allclean: clean
	rm -f $(BIN) 
