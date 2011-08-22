OBJ_DIR = ../../obj

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES = BandedSmithWaterman.cpp \
		SmithWatermanGotoh.cpp

OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES) )
INCLUDE = -I..

# ----------------
# compiler options
# ----------------

PROGRAM=smithwaterman

all: $(OBJECTS)

.PHONY: all

$(OBJECTS): $(SOURCES)
	@echo "  * compliing" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp $(INCLUDE) -o $@

