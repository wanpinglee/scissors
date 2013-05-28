# define our object and binary directories
OBJ_DIR = ../obj
BIN_DIR = ../bin

export CXX = g++
export CC = gcc
export MAKE = make
#export CFLAGS = -Wall -Werror -std=gnu99 -O3
#export CXXFLAGS = -Wall -O3
export CFLAGS = -Wall -Werror -std=gnu99 -O2 -g # enables verbose debugging
export CXXFLAGS = -Wall -Werror -O2 -g # enables verbose debugging

# define our sources
SOURCES = main.cpp
OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES) )
LIBS = -lz -lpthread
LIBS_DIR = -Loutsources/samtools
INCLUDE = -I.

EXECUTABLE = scissors

# all
all: $(EXECUTABLE)

# 
api: $(OBJECTS)
	ar rvs api/libscissors.a $(OBJ_DIR)/*.o

# test
test: $(OBJECTS)
	@cd googletest && $(MAKE) --no-print-directory

# clean
clean:
	@echo "Cleaning up......"
	@rm -rf $(OBJ_DIR) 
	@rm -rf $(BIN_DIR)
	@cd outsources/samtools && $(MAKE) clean --no-print-directory
	@rm api/libscissors.a

# define our source subdirectories
SUB_DIRS = dataStructures \
		utilities \
		outsources

SAMDIR = outsources/samtools
SAMOBJS = $(SAMDIR)/bam.o \
	  $(SAMDIR)/bgzf.o \
	  $(SAMDIR)/bam_index.o \
	  $(SAMDIR)/bam_aux.o \
	  $(SAMDIR)/sam_header.o \
	  $(SAMDIR)/bam_import.o \
	  $(SAMDIR)/kstring.o \
	  $(SAMDIR)/knetfile.o

$(EXECUTABLE): $(OBJECTS) | $(BIN_DIR)
	@echo "  * linking $(EXECUTABLE)"
	@$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/$@ $(OBJ_DIR)/*.o $(SAMOBJS) $(LIBS_DIR) $(LIBS) 

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

.PHONY: all, api, test, clean
