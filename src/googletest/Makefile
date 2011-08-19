OBJ_DIR = ../../obj
BIN_DIR = ../../bin

# define google test sources, objects, and libraries
GINCLUDES = google_include
GTEST_DIR = google_src
GTEST_HEADERS = $(GINCLUDES)/gtest/*.h \
		$(GINCLUDES)/gtest/internal/*.h

GTEST_SRCS_ = $(GTEST_DIR)/*.cc \
		$(GTEST_DIR)/*.h \
		$(GTEST_HEADERS)

GTEST_OBJS = $(OBJ_DIR)/gtest-all.o \
		$(OBJ_DIR)/gtest_main.o

GTEST_BUILD = $(GTEST_OBJS) $(OBJ_DIR)/gtest_main.a
# end //define google test sources, objects, and libraries

# ==========================
# define my own test sources
# ==========================
TESTS_DIR = test_src
TESTS_SOURCES = bam_utilities_GetPackedCigar_test.cpp
TESTS_OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(TESTS_SOURCES) )
TESTS_INCLUDES = ..
TEST_EXE = aligner_test
# end //define my own test sources

# ========================
# define my target sources
# ========================
TARGET_SOURCES = bam_utilities.cpp
TARGET_OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(TARGET_SOURCES) )
# end // define my target sources

all: $(TEST_EXE)

$(TEST_EXE): $(GTEST_BUILD) $(TESTS_OBJECTS) | $(BIN_DIR)
	@echo "  * linking" $(TEST_EXE)
	@$(CXX) -I$(GINCLUDES) -I$(GTEST_DIR) $^ $(TARGET_OBJECTS) -lpthread -o $(BIN_DIR)/$@


# Builds gtest library
$(GTEST_OBJS): $(GTEST_SRCS_)
	@echo "  * compiling" $(GTEST_DIR)/$(*F).cc
	@$(CXX) $(CXXFLAGS) -I$(GINCLUDES) -c $(GTEST_DIR)/$(*F).cc -o $@

$(OBJ_DIR)/gtest_main.a: $(GTEST_OBJS)
	@echo "  * creating gtest_main.a"
	@$(AR) -r $@ $^
# end // Builds gtest library


# Builds my own test sources
$(TESTS_OBJECTS): $(TESTS_DIR)/$(TESTS_SOURCES)
	@echo "  * compiling" $(TESTS_DIR)/$(*F).cpp
	@$(CXX) $(CXXFLAGS) -I$(TESTS_INCLUDES) -I$(GINCLUDES) -c $(TESTS_DIR)/$(*F).cpp -o $@
# end // Builds my own test sources


$(BIN_DIR):
	@mkdir $(BIN_DIR)


