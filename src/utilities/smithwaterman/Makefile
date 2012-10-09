OBJ_DIR = ../../../obj

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES = BandedSmithWaterman.cpp \
		SmithWatermanGotoh.cpp \
		ssw_cpp.cpp 

OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES) )
INCLUDE = -I../..

CSOURCES = ssw.c
COBJECTS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(CSOURCES) )

# ----------------
# compiler options
# ----------------

all: $(OBJECTS) $(COBJECTS)

.PHONY: all

$(OBJECTS): $(SOURCES)
	@echo "  * compliing" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp $(INCLUDE) -o $@

$(COBJECTS): $(CSOURCES)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c $(INCLUDE) -o $@
