CXX = g++

CXXFLAGS = -std=c++11 -fopenmp

SOURCES = $(wildcard *.cpp) $(wildcard **/*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

TARGET = main

$(TARGET) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm **/*.o *.o main