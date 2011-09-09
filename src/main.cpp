#include <stdlib.h>

#include <iostream>
#include <string>

extern "C" {
#include "dataStructures/SR_QueryRegion.h"
#include "hashTable/common/SR_Types.h"
#include "hashTable/reader/SR_HashRegionTable.h"
#include "hashTable/reader/SR_InHashTable.h"
#include "hashTable/reader/SR_Reference.h"
#include "samtools/bam.h"
#include "utilities/SR_BamInStream.h"
}

#include "dataStructures/anchor_region.h"
#include "dataStructures/search_region_type.h"
#include "dataStructures/technology.h"
#include "utilities/bam_utilities.h"
#include "utilities/parameter_parser.h"

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
  SR_QueryRegion*  query_region;
  SR_BamHeader*    bam_header;
  SR_Reference*    reference;
  SR_RefHeader*    reference_header;
  SR_InHashTable*  hash_table;
  AnchorRegion     anchor_region;
  SearchRegionType search_region_type;
  RegionType       region_type;
  BamReference     bam_reference;
};


// Prototype of functions
void Deconstruct(MainFiles* files, MainVars* vars);
void PrepareFiles(const ParameterParser& parameter_parser, MainFiles* files);
void PrepareVariablesOrDie(
    const ParameterParser& parameter_parser, 
    const MainFiles& files, MainVars* vars);
void CheckFileOrDie(
    const ParameterParser& parameter_parser,
    const MainFiles& files);
void ResetHeader( bam_header_t* const bam_header );


int main ( int argc, char** argv ) {
	
  // Parse the arguments and store them
  // The program will exit(1) with printing error message 
  //     if any errors or missing required parameters are found
  ParameterParser const parameter_parser(argc, argv);


  // =================
  // Files Preparation
  // =================
  MainFiles files;
  PrepareFiles(parameter_parser, &files);
  CheckFileOrDie(parameter_parser, files);


  // =====================
  // Variables Preparation
  // =====================
  MainVars vars;
  PrepareVariablesOrDie(parameter_parser, files, &vars);

  // Write bam header
  ResetHeader( vars.bam_header->pOrigHeader );
  bam_header_write( files.bam_writer, vars.bam_header->pOrigHeader );


  // Load first reference and its hash table
  SR_ReferenceRead(vars.reference, files.ref_reader);
  SR_InHashTableRead(vars.hash_table, files.hash_reader);


  // =========
  // Algorithm
  // =========

  while(SR_BamInStreamGetPair( &(vars.query_region->pAnchor), &(vars.query_region->pOrphan), files.bam_reader ) == SR_OK) {
    uint32_t* cigar = bam1_cigar(vars.query_region->pAnchor);
    bool is_new_region = vars.anchor_region.IsNewRegion(cigar, 
                         vars.query_region->pAnchor->core.n_cigar, vars.query_region->pAnchor->core.pos);
		
    if (is_new_region) {
      vars.search_region_type.ResetRegionTypeList();
    }else{
      vars.search_region_type.RewindRegionTypeList();
    }

    const bool is_anchor_forward = bam1_strand(vars.query_region->pAnchor);
    while (vars.search_region_type.GetNextRegionType(is_anchor_forward, &vars.region_type))
      printf("-");
    cout << "Got a pair of alignments" << endl;
		
    //bam_write1( files.bam_writer, vars.query_region->pAnchor );
    //bam_write1( files.bam_writer, vars.query_region->pOrphan );
  }

  // free memory and close files
  Deconstruct(&files, &vars);

  cout << "Program done." << endl;

  return 0;

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

}

void PrepareFiles(const ParameterParser& parameter_parser, MainFiles* files) {
  // Initialize bam input reader
  // The program will be terminated with printing error message
  //     if the given bam cannot be opened.
  files->bam_reader = SR_BamInStreamAlloc(
      parameter_parser.input_bam.c_str(), 
      parameter_parser.fragment_length * parameter_parser.mate_window_size, 
      parameter_parser.allowed_clip );
  // Initialize bam output writer
  files->bam_writer = bam_open( parameter_parser.output_bam.c_str(), "w" );

  // Initialize reference input reader
  files->ref_reader  = fopen( parameter_parser.reference_filename.c_str(), "rb");
  files->hash_reader = fopen( parameter_parser.hash_filename.c_str(), "rb");

}

void PrepareVariablesOrDie(const ParameterParser& parameter_parser, 
    const MainFiles& files,
    MainVars* vars) {
  // Load bam header	
  vars->bam_header = SR_BamHeaderAlloc();
  vars->bam_header = SR_BamInStreamLoadHeader(files.bam_reader);
  if( !parameter_parser.is_input_sorted && !BamUtilities::IsFileSorted( vars->bam_header->pOrigHeader ) ) {
    // The input bam is unsorted, exit
    cout << "ERROR: The input bam seems unsorted. "
         << "Please use bamtools sort to sort the bam" << endl
	 << "       or type -s to ignore this checker." << endl;
	 exit(1);
  }

  // bam reference
  vars->bam_reference.num_reference = 
    SR_BamHeaderGetRefNum(vars->bam_header);
  vars->bam_reference.reference_names = 
    SR_BamHeaderGetRefNames(vars->bam_header);
  vars->bam_reference.reference_lengths = 
    SR_BamHeaderGetRefLens(vars->bam_header);

  // Note: bam records are in SR_QueryRegion structure
  vars->query_region = SR_QueryRegionAlloc();

  // Load header of reference file
  vars->reference_header  = SR_RefHeaderAlloc();
  int64_t reference_seal  = SR_RefHeaderRead(vars->reference_header, files.ref_reader);

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
}

void CheckFileOrDie( 
		const ParameterParser& parameter_parser,
		const MainFiles& files){

	bool error_found = false;

	if ( files.bam_writer == NULL ) {
		cout << "ERROR: Cannot open " << parameter_parser.output_bam << " for writing." << endl
		     << "       Please check -o option." << endl;
		error_found = true;
	}

	if ( files.ref_reader == NULL ) {
		cout << "ERROR: Cannot open " << parameter_parser.reference_filename << " for reading." << endl
		     << "       Please check -r option." << endl;
		error_found = true;
	}

	if ( files.hash_reader == NULL ) {
		cout << "ERROR: Cannot open " << parameter_parser.hash_filename << " for reading." << endl
		     << "       Please check -r option." << endl;
		error_found = true;
	}

	if ( error_found )
		exit(1);

}

void ResetHeader( bam_header_t* const bam_header ){

	// Reset header line to "@HD\tVN:1.0\tSO:unsorted"
	//BamUtilities::ResetHeaderLineText( bam_header, "@HD\tVN:1.0\tSO:unsorted" );
	// Replace SO:coordinate or SO:queryname by SO:unsorted
	BamUtilities::ReplaceHeaderSoText( bam_header );
}

