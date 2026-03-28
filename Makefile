CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic -Isrc
LDFLAGS  = -lGL -lGLU -lglut -lm

SRC = src/main.cpp \
      src/scene.cpp \
      src/hemicube.cpp \
      src/radiosity.cpp \
      src/render.cpp

OBJ    = $(SRC:.cpp=.o)
TARGET = radiosity

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(TARGET)
