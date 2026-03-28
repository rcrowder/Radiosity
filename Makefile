CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic -Isrc
LDFLAGS  = -lGL -lGLU -lglut -lm

SRC = src/main.cpp \
      src/scene.cpp \
      src/hemicube.cpp \
      src/radiosity.cpp \
      src/render.cpp \
      src/discmesh.cpp

OBJDIR = obj
OBJ    = $(patsubst src/%.cpp,$(OBJDIR)/%.o,$(SRC))
TARGET = radiosity

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET)
