OBJ_DIR = ../../obj

# C++
SOURCES = search_region_type.cpp \
		anchor_region.cpp 
OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
INCLUDE = -I..

# C
CSOURCES = SR_QueryRegion.c
COBJECTS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(CSOURCES))
CINCLUDE = -I..

all: $(OBJECTS) $(COBJECTS)

$(OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp $(INCLUDE) -o $@

$(COBJECTS): $(CSOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c $(CINCLUDE) -o $@

.PHONY: all

