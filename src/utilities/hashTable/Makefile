OBJ_DIR = ../../../obj

# C++
SOURCES = special_hasher.cpp \
		reference_hasher.cpp 

OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
INCLUDE = -I../..

CSOURCES = SR_InHashTable.c\
		SR_OutHashTable.c \
		SR_Reference.c \
		SR_HashRegionTable.c \
		ConvertHashTableOutToIn.c

COBJECTS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(CSOURCES) )
CINCLUDE = -I../..

all: $(OBJECTS) $(COBJECTS)

$(OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp $(INCLUDE) -o $@

$(COBJECTS): $(CSOURCES)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c $(CINCLUDE) -o $@

.PHONY: all
