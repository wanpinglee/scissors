#include <stdlib.h>

#include <iostream>
#include <string>

extern "C" {
#include "outsources/samtools/bam.h"
#include "utilities/bam/SR_BamInStream.h"
#include "utilities/bam/SR_BamPairAux.h"
#include "utilities/common/SR_Types.h"
}

#include "utilities/bam/bam_reference.h"
#include "utilities/bam/bam_utilities.h"
//#include "utilities/bam/bam_writer.h"
#include "utilities/miscellaneous/parameter_parser.h"
#include "utilities/miscellaneous/thread.h"

using std::string;
using std::cout;
using std::endl;

struct MainFiles {
  SR_BamInStream* bam_reader;  // bam reader
  bamFile         bam_writer;  // bam writer
  FILE*           ref_reader;  // reference reader
  FILE*           hash_reader; // hash table reader
};

struct MainVars{
  // for bam
  SR_BamHeader*    bam_header;
  SR_StreamMode    bam_record_mode;
};


// Prototype of functions
void Deconstruct(MainFiles* files, MainVars* vars);
void InitFiles(const Parameters& parameters, MainFiles* files);
void InitVariablesOrDie(const Parameters& parameter_parser, 
                        const MainFiles& files, 
			MainVars* vars);
void CheckFileOrDie(const Parameters& parameters,
                    const MainFiles& files);
void ResetHeader(bam_header_t* const bam_header);


int main ( int argc, char** argv ) {
	
  // Parse the arguments and store them
  // The program will exit(1) with printing error message 
  //     if any errors or missing required parameters are found
  Parameters parameters;
  ParseArgumentsOrDie(argc, argv, &parameters);

  // =================
  // Files Preparation
  // =================
  MainFiles files;
  InitFiles(parameters, &files);
  CheckFileOrDie(parameters, files);


  // =====================
  // Variables Preparation
  // =====================
  MainVars vars;
  InitVariablesOrDie(parameters, files, &vars);

  // Write bam header
  ResetHeader(vars.bam_header->pOrigHeader);
  bam_header_write( files.bam_writer, vars.bam_header->pOrigHeader );

  // =========
  // Algorithm
  // =========
  BamReference bam_reference;
  bam_reference.Init(*vars.bam_header);
  Thread thread(&bam_reference,
		parameters.allowed_clip,
		parameters.processors,
		files.ref_reader,
		files.hash_reader,
		files.bam_reader,
		&files.bam_writer);
  bool thread_status = thread.Start();
  if (!thread_status)
    cout << "threads fail" << endl;
  //StartThreadOrDie(parameters.processors, files.bam_reader);

  // free memory and close files
  Deconstruct(&files, &vars);

  cout << "Program done." << endl;

  return 0;

}

/*
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
*/

void Deconstruct(MainFiles* files, MainVars* vars) {
  // close files
  SR_BamInStreamFree(files->bam_reader);
  bam_close(files->bam_writer);
  //files->bam_writer.Close();
  fclose(files->ref_reader);
  fclose(files->hash_reader);

  // free variables
  SR_BamHeaderFree(vars->bam_header);
}

void InitFiles(const Parameters& parameters, MainFiles* files) {
  // set the stream mode to "UO" (unique orphan)
  SR_StreamMode streamMode;
  SR_SetStreamMode(&streamMode, SR_CommonFilter, NULL, SR_NO_SPECIAL_CONTROL);
  // Initialize bam input reader
  // The program will be terminated with printing error message
  //     if the given bam cannot be opened.
  files->bam_reader = SR_BamInStreamAlloc(
      parameters.input_bam.c_str(), 
      parameters.fragment_length * parameters.mate_window_size,
      parameters.processors,  // number of processors
      2, // the number of alignments can be stored in each chunk of the memory pool
      2, // number of alignments should be cached before report
      &streamMode);

  // Initialize bam output writer
  files->bam_writer = bam_open( parameters.output_bam.c_str(), "w" );
  //files->bam_writer.Open(parameters.output_bam);

  // Initialize reference input reader
  files->ref_reader  = fopen( parameters.reference_filename.c_str(), "rb");
  files->hash_reader = fopen( parameters.hash_filename.c_str(), "rb");

}

void IsInputBamSortedOrDie(const Parameters& parameters,
                           const SR_BamHeader& bam_header) {
  if (!parameters.is_input_sorted &&
      !BamUtilities::IsFileSorted(bam_header.pOrigHeader)) {
    // The input bam is unsorted, exit
    cout << "ERROR: The input bam seems unsorted. "
         << "Please use bamtools sort to sort the bam" << endl
	 << "       or type -s to ignore this checker." << endl;
	 exit(1);
  }
}

void InitVariablesOrDie(const Parameters& parameters, 
                        const MainFiles& files,
                        MainVars* vars) {
  // Load bam header	
  vars->bam_header = SR_BamHeaderAlloc();
  vars->bam_header = SR_BamInStreamLoadHeader(files.bam_reader);

  IsInputBamSortedOrDie(parameters, *(vars->bam_header));

  //vars->search_window.fragLen    = 1000;
  //vars->search_window.closeRange = 2000;
  //vars->search_window.farRange   = 100000;
}

void CheckFileOrDie(const Parameters& parameters,
                    const MainFiles& files){

	bool error_found = false;

	if (files.bam_writer == NULL) {
		cout << "ERROR: Cannot open " 
		     << parameters.output_bam 
		     << " for writing." 
		     << endl
		     << "       Please check -o option." 
		     << endl;
		error_found = true;
	}

	if (files.ref_reader == NULL) {
		cout << "ERROR: Cannot open " 
		     << parameters.reference_filename 
		     << " for reading." 
		     << endl
		     << "       Please check -r option." << endl;
		error_found = true;
	}

	if (files.hash_reader == NULL) {
		cout << "ERROR: Cannot open " 
		     << parameters.hash_filename 
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

