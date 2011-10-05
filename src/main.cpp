#include <stdlib.h>

#include <iostream>
#include <string>

extern "C" {
#include "dataStructures/SR_QueryRegion.h"
#include "outsources/samtools/bam.h"
#include "utilities/bam/SR_BamHeader.h"
#include "utilities/bam/SR_BamInStream.h"
#include "utilities/bam/SR_BamPairAux.h"
#include "utilities/common/SR_Types.h"
#include "utilities/hashTable/SR_HashRegionTable.h"
#include "utilities/hashTable/SR_InHashTable.h"
#include "utilities/hashTable/SR_Reference.h"
}

#include "dataStructures/anchor_region.h"
#include "dataStructures/search_region_type.h"
#include "dataStructures/technology.h"
#include "utilities/bam/bam_utilities.h"
#include "utilities/miscellaneous/hash_region_collection.h"
#include "utilities/miscellaneous/parameter_parser.h"

using std::string;
using std::cout;
using std::endl;

struct MainFiles {
  SR_BamInStream* bam_reader;  // bam reader
  bamFile         bam_writer;  // bam writer
  FILE*           ref_reader;  // reference reader
  FILE*           hash_reader; // hash table reader
};

struct BamReference {
  int             num_reference;
  const char**    reference_names;
  const uint32_t* reference_lengths;
};

struct MainVars{
  // for bam
  SR_BamHeader*    bam_header;
  SR_StreamMode    bam_record_mode;
  BamReference     bam_reference;

  SR_QueryRegion*  query_region;
  SR_Reference*    reference;
  SR_RefHeader*    reference_header;
  SR_InHashTable*  hash_table;

  AnchorRegion     anchor_region;
  SearchRegionType search_region_type;
  HashRegionTable* hash_region_table;
  SR_SearchArgs    search_window;
  SearchRegionType::RegionType       region_type;

  HashRegionCollection hr_collection;

  float soft_clip_tolerance;
};


// Prototype of functions
void Deconstruct(MainFiles* files, MainVars* vars);
void InitFiles(const ParameterParser& parameter_parser, MainFiles* files);
void InitVariablesOrDie(const ParameterParser& parameter_parser, 
                        const MainFiles& files, 
			MainVars* vars);
void LoadReferenceOrDie(const bam1_t& anchor,
                        MainFiles* files, MainVars* vars);
void CheckFileOrDie(const ParameterParser& parameter_parser,
                    const MainFiles& files);
void ResetHeader(bam_header_t* const bam_header);
void LoadRegionType(const bam1_t& anchor,
                    AnchorRegion& anchor_region,
                    SearchRegionType& search_region_type);
void SetTargetSequence(const SearchRegionType::RegionType& region_type,
                       SR_QueryRegion* query_region);
void AlignCandidate(const SR_BamListIter& al_ite,
                    MainFiles* files,
                    MainVars* vars); 


int main ( int argc, char** argv ) {
	
  // Parse the arguments and store them
  // The program will exit(1) with printing error message 
  //     if any errors or missing required parameters are found
  Parameters param;
  ParseArgumentsOrDie(argc, argv, &param);

  // =================
  // Files Preparation
  // =================
  MainFiles files;
  InitFiles(parameter_parser, &files);
  CheckFileOrDie(parameter_parser, files);


  // =====================
  // Variables Preparation
  // =====================
  MainVars vars;
  InitVariablesOrDie(parameter_parser, files, &vars);

  // Write bam header
  ResetHeader( vars.bam_header->pOrigHeader );
  bam_header_write( files.bam_writer, vars.bam_header->pOrigHeader );

  // =========
  // Algorithm
  // =========
  SR_Status bam_status;
  const int threads = 2;
  vector<SR_BamListIter> iters;
  iters.resize(threads);

  for (int i = 0; i < threads; ++i) {
    // try to load the first chunk of alignments
    bam_status = SR_LoadUniquOrphanPairs(files.bam_reader, i, vars.soft_clip_tolerance);
  }

  if (bam_status != SR_ERR) {
    for (int threadId = 0; threadId < threads; ++threadId) {
      iters[threadId] = SR_BamInStreamGetIter(files.bam_reader, threadId);
    }
  } else {
    cout << "Not thing is loaded." << endl;
  } // end if-else

		
  // free memory and close files
  Deconstruct(&files, &vars);

  cout << "Program done." << endl;

  return 0;

}

void LoadRegionType(const bam1_t& anchor, 
    AnchorRegion& anchor_region,
    SearchRegionType& search_region_type) {
  uint32_t* cigar = bam1_cigar(&anchor);
  bool is_new_region = anchor_region.IsNewRegion(cigar,
                       anchor.core.n_cigar, 
		       anchor.core.pos);
  if (is_new_region)
    search_region_type.ResetRegionTypeList();
  else
    search_region_type.RewindRegionTypeList();

}

void Deconstruct(MainFiles* files, MainVars* vars) {
  // close files
  SR_BamInStreamFree(files->bam_reader);
  bam_close(files->bam_writer);
  fclose(files->ref_reader);
  fclose(files->hash_reader);

  // free variables
  SR_QueryRegionFree(vars->query_region);
  SR_BamHeaderFree(vars->bam_header);
  SR_ReferenceFree(vars->reference);
  SR_InHashTableFree(vars->hash_table);
  HashRegionTableFree(vars->hash_region_table);

}

void InitFiles(const ParameterParser& parameter_parser, MainFiles* files) {
  // set the stream mode to "UO" (unique orphan)
  SR_StreamMode streamMode;
  SR_SetStreamMode(&streamMode, SR_UO_STREAM);
  // Initialize bam input reader
  // The program will be terminated with printing error message
  //     if the given bam cannot be opened.
  files->bam_reader = SR_BamInStreamAlloc(
      parameter_parser.input_bam.c_str(), 
      parameter_parser.fragment_length * parameter_parser.mate_window_size,
      2,  // number of threads
      20, // the number of alignments can be stored in each chunk of the memory pool
      200, // number of alignments should be cached before report
      &streamMode);

  // Initialize bam output writer
  files->bam_writer = bam_open( parameter_parser.output_bam.c_str(), "w" );

  // Initialize reference input reader
  files->ref_reader  = fopen( parameter_parser.reference_filename.c_str(), "rb");
  files->hash_reader = fopen( parameter_parser.hash_filename.c_str(), "rb");

}

void InitBamReference(const SR_BamHeader& bam_header,
                      BamReference* bam_reference) {
  bam_reference->num_reference     = SR_BamHeaderGetRefNum(&bam_header);
  bam_reference->reference_names   = SR_BamHeaderGetRefNames(&bam_header);
  bam_reference->reference_lengths = SR_BamHeaderGetRefLens(&bam_header);
}

void IsInputBamSortedOrDie(const ParameterParser& parameter_parser,
                           const SR_BamHeader& bam_header) {
  if (!parameter_parser.is_input_sorted &&
      !BamUtilities::IsFileSorted(bam_header.pOrigHeader)) {
    // The input bam is unsorted, exit
    cout << "ERROR: The input bam seems unsorted. "
         << "Please use bamtools sort to sort the bam" << endl
	 << "       or type -s to ignore this checker." << endl;
	 exit(1);
  }
}

void InitVariablesOrDie(const ParameterParser& parameter_parser, 
                        const MainFiles& files,
                        MainVars* vars) {
  // Load bam header	
  vars->bam_header = SR_BamHeaderAlloc();
  vars->bam_header = SR_BamInStreamLoadHeader(files.bam_reader);

  IsInputBamSortedOrDie(parameter_parser, *(vars->bam_header));

  // bam reference
  InitBamReference(*vars->bam_header, &vars->bam_reference);

  // Note: bam records are in SR_QueryRegion structure
  vars->query_region = SR_QueryRegionAlloc();

  // Load header of reference file
  vars->reference_header  = SR_RefHeaderAlloc();
  
  int64_t reference_seal  = 
    SR_RefHeaderRead(vars->reference_header, files.ref_reader);


  // Load header of hash file
  unsigned char hash_size = 0;
  int64_t hash_seal = SR_InHashTableReadStart(&hash_size, files.hash_reader);

  if (reference_seal != hash_seal) {
    printf("ERROR: The reference file is not compatible with the hash table file.\n");
    exit(1);
  }

  // init reference and hash table
  vars->reference  = SR_ReferenceAlloc();
  vars->hash_table = SR_InHashTableAlloc(hash_size);

  vars->hash_region_table = HashRegionTableAlloc();

  vars->search_window.fragLen    = 1000;
  vars->search_window.closeRange = 2000;
  vars->search_window.farRange   = 100000;
  vars->soft_clip_tolerance      = 0.2;
}


void LoadReferenceOrDie(const bam1_t& anchor,
                        MainFiles* files, MainVars* vars) {
  // Unknown reference id
  if (anchor.core.tid > vars->bam_reference.num_reference) {
    printf("ERROR: The reference id of the anchor, %s, is invalid.\n", bam1_qname(&anchor));
    exit(1);
  }

  const char* bam_ptr = vars->bam_reference.reference_names[anchor.core.tid];
  const char* ref_ptr = SR_RefHeaderGetName(vars->reference_header, vars->reference->id);
  
  if (strcmp(bam_ptr, ref_ptr) != 0) {  // Loading reference and hash table is necessary
      int32_t ref_id = SR_RefHeaderGetRefID(vars->reference_header, ref_ptr);
      
      if (ref_id < 0) {  // Cannot find the reference
        printf("ERROR: The reference, %s, is not found in reference file.\n", ref_ptr);
	exit(1);
      } else {
        SR_ReferenceJump(files->ref_reader, vars->reference_header, ref_id);
        SR_InHashTableJump(files->hash_reader, vars->reference_header, ref_id);
        SR_ReferenceRead(vars->reference, files->ref_reader);
        SR_InHashTableRead(vars->hash_table, files->hash_reader);
      }
  }

}

void CheckFileOrDie(const ParameterParser& parameter_parser,
                    const MainFiles& files){

	bool error_found = false;

	if (files.bam_writer == NULL) {
		cout << "ERROR: Cannot open " 
		     << parameter_parser.output_bam 
		     << " for writing." 
		     << endl
		     << "       Please check -o option." 
		     << endl;
		error_found = true;
	}

	if (files.ref_reader == NULL) {
		cout << "ERROR: Cannot open " 
		     << parameter_parser.reference_filename 
		     << " for reading." 
		     << endl
		     << "       Please check -r option." << endl;
		error_found = true;
	}

	if (files.hash_reader == NULL) {
		cout << "ERROR: Cannot open " 
		     << parameter_parser.hash_filename 
		     << " for reading." 
		     << endl
		     << "       Please check -r option." 
		     << endl;
		error_found = true;
	}

	if (error_found)
		exit(1);

}

void ResetHeader( bam_header_t* const bam_header ){

	// Reset header line to "@HD\tVN:1.0\tSO:unsorted"
	//BamUtilities::ResetHeaderLineText( bam_header, "@HD\tVN:1.0\tSO:unsorted" );
	// Replace SO:coordinate or SO:queryname by SO:unsorted
	BamUtilities::ReplaceHeaderSoText( bam_header );
}

