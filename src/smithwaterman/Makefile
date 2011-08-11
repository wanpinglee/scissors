OBJ_DIR = ../../obj

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES = BandedSmithWaterman.cpp \
		SmithWatermanGotoh.cpp

OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS = $(patsubst %, $(OBJ_DIR)/%, $(OBJECTS) )

# ----------------
# compiler options
# ----------------

PROGRAM=smithwaterman

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compliing" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp -o $@

