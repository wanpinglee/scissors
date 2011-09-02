OBJ_DIR = ../../obj
BIN_DIR = ../../bin

TESTS_DIR = test_src
# *************************************
# * PLEASE ADD GOOGLE TEST FILES HERE *
# * THE FILES SHOULD BE IN TESTS_DIR  *
# *************************************
TESTS_SOURCES_ = bam_utilities_GetPackedCigar_test.cpp \
		anchor_region_test.cpp \
		search_region_type_test.cpp

TARGET_SOURCES = bam_utilities.cpp \
			anchor_region.cpp \
			search_region_type.cpp

# *********************************************
# YOU MAY NOT NEED TO CHANGE THE FOLLOWING CODE
# *********************************************

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
TESTS_SOURCES = $(patsubst %, $(TESTS_DIR)/%, $(TESTS_SOURCES_) )
TESTS_OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(TESTS_SOURCES_) )
TESTS_INCLUDES = ..
TEST_EXE = aligner_test
# end //define my own test sources

# ========================
# define my target sources
# ========================
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
$(TESTS_OBJECTS): $(TESTS_SOURCES)
	@echo "  * compiling" $(TESTS_DIR)/$(*F).cpp
	@$(CXX) $(CXXFLAGS) -I$(TESTS_INCLUDES) -I$(GINCLUDES) -c $(TESTS_DIR)/$(*F).cpp -o $@
# end // Builds my own test sources


$(BIN_DIR):
	@mkdir $(BIN_DIR)


