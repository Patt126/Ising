CXX = mpicxx
CXXFLAGS = -I../muparser/include -L../muparser/build -lmuparser -Wall -fopenmp
SRC = main.cpp MontIntegration.hpp Shape.hpp HyperRectangle.hpp HyperSphere.hpp MontIntegration.hpp
OBJ = $(SRC:.cpp=.o)
EXEC = main

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(OBJ) -o $(EXEC) $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(OBJ) $(EXEC)
