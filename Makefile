CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic -Isrc
LDFLAGS  = -lGL -lGLU -lglut -lm

SRC = src/main.cpp \
      src/scene.cpp \
      src/hemicube.cpp \
      src/radiosity.cpp \
      src/render.cpp \
      src/discmesh.cpp \
      src/teapot.cpp

OBJDIR = obj
OBJ    = $(patsubst src/%.cpp,$(OBJDIR)/%.o,$(SRC))
TARGET = radiosity

TEST_TARGET = tests/test_discmesh
TEST_SRC    = tests/test_discmesh.cpp src/scene.cpp src/discmesh.cpp

.PHONY: all clean test

all: $(TARGET)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

$(TEST_TARGET): $(TEST_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET)
