OBJ_DIR = ../../obj

SOURCES = ParameterParser.cpp
OBJECTS = $(SOURCES:.cpp=.o)
BUILT_OBJECTS = $(patsubst %, $(OBJ_DIR)/%, $(OBJECTS) )

all: $(BUILT_OBJECTS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compliing" $(*F).cpp
	@$(CXX) -c $(CXXFLAGS) $(*F).cpp -o $@

.PHONY: all

