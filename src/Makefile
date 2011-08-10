# define our object and binary directories
OBJ_DIR = ../obj
BIN_DIR = ../bin

export CXX = g++
export CC = gcc
export CFLAGS = -Wall -O3 -std=c99
export CXXFLAGS = -Wall -O3
export MAKE = make

# define our sources
SOURCES = main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
BUILT_OBJECTS = $(OBJECTS:%=$(OBJ_DIR)/%)

# define our source subdirectories
SUB_DIRS = hasher utilities
ext_objects = ParameterParser.o
EXT_OBJECTS = $(ext_objects:%=$(OBJ_DIR)/%)

PROGRAM = aligner

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(EXT_OBJECTS) $(BUILT_OBJECTS)
	@test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ 

$(BUILT_OBJECTS): $(SOURCES)
	@test -d $(OBJ_DIR) || mkdir $(OBJ_DIR)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS)

$(EXT_OBJECTS):
	@echo "Building......"
	@for dir in $(SUB_DIRS); do \
		$(MAKE) --no-print-directory --directory=$$dir; \
	done

clean:
	@echo "Cleaning......"
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
.PHONY: clean
