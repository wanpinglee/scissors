OBJ_DIR = ../../obj
BIN_DIR = ../../bin

TESTS_DIR = test_src
# *************************************
# * PLEASE ADD GOOGLE TEST FILES HERE *
# * THE FILES SHOULD BE IN TESTS_DIR  *
# *************************************
TESTS_SOURCES_ = bam_utilities_GetPackedCigar_test.cpp \
		anchor_region_test.cpp \
		search_region_type_test.cpp \
		aligner_api_test.cpp \
		alignment_filter_test.cpp

TARGET_OBJECTS_ = bam_utilities.o \
			anchor_region.o \
			search_region_type.o \
			special_hasher.o \
			reference_hasher.o \
			ConvertHashTableOutToIn.o \
			aligner.o \
			BandedSmithWaterman.o \
			hashes_collection.o \
			optional_tag.o \
			alignment_filter.o \
			SR_QueryRegion.o \
			SR_HashRegionTable.o \
			SR_BamPairAux.o \
			SR_BamInStream.o \
			SR_BamMemPool.o \
			SR_FragLenDstrb.o \
			SR_BamHeader.o


REQUIRED_OBJS_ = SR_Reference.o \
		SR_OutHashTable.o \
		SR_InHashTable.o \
		SR_Error.o

SAMDIR = ../outsources/samtools
SAMOBJS_ = misc/md5.o \
		bam.o \
		bgzf.o \
		bam_index.o \
		bam_aux.o \
		sam_header.o \
		bam_import.o \
		kstring.o \
		knetfile.o


AUXDIR = util
AUX_SOURCES_ = fasta_reader.cpp

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
AUX_SOURCES = $(patsubst %.cpp, $(AUXDIR)/%.cpp, $(AUX_SOURCES_) )
TESTS_SOURCES = $(patsubst %, $(TESTS_DIR)/%, $(TESTS_SOURCES_) )
TESTS_OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(TESTS_SOURCES_) )
REQUIRED_OBJS = $(patsubst %.o, $(OBJ_DIR)/%.o, $(REQUIRED_OBJS_) )
SAMOBJS = $(patsubst %.o, $(SAMDIR)/%.o, $(SAMOBJS_))
AUX_OBJS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(AUX_SOURCES_) )
TESTS_INCLUDES = ..
TEST_EXE = aligner_test
# end //define my own test sources

# ========================
# define my target sources
# ========================
TARGET_OBJECTS = $(patsubst %.o, $(OBJ_DIR)/%.o, $(TARGET_OBJECTS_) )
# end // define my target sources

all: $(TEST_EXE)

$(TEST_EXE): $(GTEST_BUILD) $(TESTS_OBJECTS) $(AUX_OBJS) | $(BIN_DIR)
	@echo "  * linking" $(TEST_EXE)
	@$(CXX) -I$(GINCLUDES) -I$(GTEST_DIR) $^ $(TARGET_OBJECTS) $(REQUIRED_OBJS) $(SAMOBJS) -lz -lpthread -o $(BIN_DIR)/$@


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

$(AUX_OBJS): $(AUX_SOURCES)
	@echo "  * compiling" $(AUX_SOURCES)
	@$(CXX) $(CXXFLAGS) -c $(AUX_SOURCES) -o $@

$(BIN_DIR):
	@mkdir $(BIN_DIR)


