OBJ_DIR = ../../obj

SOURCES = ParameterParser.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(OBJECTS)

$(OBJECTS): $(SOURCES)
	@test -d $(OBJ_DIR) || mkdir $(OBJ_DIR)
	@echo "  * compliing" $(*F).cpp
	@$(CXX) -c -o $(OBJ_DIR)/$@ $(*F).cpp $(CXXFLAGS)

.PHONY: all

