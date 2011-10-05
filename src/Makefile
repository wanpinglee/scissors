# define our object and binary directories
OBJ_DIR = ../obj
BIN_DIR = ../bin

export CXX = g++
export CC = gcc
export MAKE = make
export CFLAGS = -Wall -std=gnu99 -g
export CXXFLAGS = -Wall -g

# define our sources
SOURCES = main.cpp
OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES) )
LIBS = -lz -lbam
LIBS_DIR = -Loutsources/samtools
INCLUDE = -I.


# define our source subdirectories
SUB_DIRS = dataStructures \
		utilities \
		outsources

EXECUTABLE = aligner

all: $(EXECUTABLE)
.PHONY: all

$(EXECUTABLE): $(OBJECTS) | $(BIN_DIR)
	@echo "  * linking $(EXECUTABLE)"
	@$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/$@ $(OBJ_DIR)/*.o $(LIBS_DIR) $(LIBS) 

test: $(OBJECTS)
	@cd googletest && $(MAKE) --no-print-directory

$(OBJECTS): $(SOURCES) | $(OBJ_DIR)
	@for dir in $(SUB_DIRS); do \
	$(MAKE) --no-print-directory --directory=$$dir; \
	done
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp $(INCLUDE) -o $@

$(OBJ_DIR):
	@mkdir $(OBJ_DIR)
$(BIN_DIR):
	@mkdir $(BIN_DIR)

clean:
	@echo "Cleaning up......"
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
	@cd outsources/samtools && $(MAKE) clean --no-print-directory
.PHONY: clean
