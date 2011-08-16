# define our object and binary directories
OBJ_DIR = ../obj
BIN_DIR = ../bin

export CXX = g++
export CC = gcc
export MAKE = make
export CFLAGS = -Wall -O3 -std=gnu99
export CXXFLAGS = -Wall -O3

# define our sources
SOURCES = main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
BUILT_OBJECTS = $(patsubst %, $(OBJ_DIR)/%, $(OBJECTS) )
LIBS = -lz


# define our source subdirectories
SUB_DIRS = hasher \
		utilities \
		smithwaterman

EXECUTABLE = aligner

all: $(EXECUTABLE)
.PHONY: all

$(EXECUTABLE): $(BUILT_OBJECTS) | $(BIN_DIR)
	@echo "  * linking $(EXECUTABLE)"
	@$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/$@ $(OBJ_DIR)/*.o $(LIBS)

$(BUILT_OBJECTS): $(SOURCES) | $(OBJ_DIR)
	@for dir in $(SUB_DIRS); do \
	$(MAKE) --no-print-directory --directory=$$dir; \
	done
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp -o $@


$(OBJ_DIR):
	@mkdir $(OBJ_DIR)
$(BIN_DIR):
	@mkdir $(BIN_DIR)

clean:
	@echo "Cleaning......"
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
.PHONY: clean
