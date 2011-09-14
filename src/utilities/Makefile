OBJ_DIR = ../../obj

# C++
SOURCES = bam_utilities.cpp \
		hash_region_collection.cpp \
		parameter_parser.cpp

OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES) )
LIBS = -lz
INCLUDE = -I..

# C
CSOURCES = md5.c \
		SR_BamInStream.c
COBJECTS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(CSOURCES) )
CINCLUDE = -I..

all: $(OBJECTS) $(COBJECTS)

$(OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp $(INCLUDE) -o $@ $(LIBS)

$(COBJECTS): $(CSOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c $(CINCLUDE) -o $@

.PHONY: all

