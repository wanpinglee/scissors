#include <stdlib.h>

#include <iostream>
#include <string>

extern "C" {
#include "outsources/samtools/bam.h"
#include "utilities/bam/SR_BamInStream.h"
#include "utilities/bam/SR_BamPairAux.h"
#include "utilities/common/SR_Types.h"
#include "utilities/miscellaneous/md5.h"
}

#include "dataStructures/target_event.h"
#include "dataStructures/target_region.h"
#include "outsources/fasta/Fasta.h"
#include "utilities/bam/bam_reference.h"
#include "utilities/bam/bam_utilities.h"
#include "utilities/miscellaneous/alignment_filter.h"
#include "utilities/miscellaneous/parameter_parser.h"
#include "utilities/miscellaneous/thread.h"

using std::string;
using std::cerr;
using std::endl;

using namespace Scissors;

struct MainFiles {
  SR_BamInStream* bam_reader;  // bam reader
  bamFile         bam_writer;  // bam writer
  FastaReference  ref_reader;
};

struct MainVars{
  // for bam
  SR_BamHeader*    bam_header;
  SR_StreamMode    bam_record_mode;
  // alignment
  AlignmentFilter  alignment_filter;
  TargetEvent      target_event;
  TargetRegion     target_region;
};


// Prototype of functions
void Deconstruct(MainFiles* files, MainVars* vars);
void InitFiles(const Parameters& parameters, MainFiles* files);
void InitVariablesOrDie(const Parameters& parameter_parser, 
                        MainFiles* files, 
			MainVars* vars);
void CheckFileOrDie(const Parameters& parameters,
                    const MainFiles& files);
void ResetSoBamHeader(bam_header_t* const bam_header);
void AppendReferenceSequence(bam_header_t* const bam_header,
                             const string& reference_filename);
void SetAlignmentFilter(const Parameters& parameters, 
                        AlignmentFilter* filter);
void SetTargetEvent(const Parameters& parameters,
                    TargetEvent* target_event);
void SetTargetRegion(const Parameters& parameters,
                     TargetRegion* target_region);


int main (int argc, char** argv) {
	
  // Parse the arguments and store them.
  // The program will exit(1) with printing error message 
  // if any errors or missing required parameters are found
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
  InitVariablesOrDie(parameters, &files, &vars);

  // Write bam header
  ResetSoBamHeader(vars.bam_header->pOrigHeader);
  if (parameters.detect_special)
    AppendReferenceSequence(vars.bam_header->pOrigHeader, parameters.input_special_fasta);
  // load the reference header
  bam_header_write(files.bam_writer, vars.bam_header->pOrigHeader);

  // =========
  // Algorithm
  // =========
  BamReference bam_reference;
  bam_reference.Init(*vars.bam_header);
  
  Thread thread(&bam_reference,
		parameters.allowed_clip,
		parameters.processors,
		parameters.fragment_length,
		parameters.technology,
		vars.target_event,
		parameters.mapping_quality_threshold,
		vars.alignment_filter,
		vars.target_region,
		parameters.input_special_fasta,
		&files.ref_reader,
		files.bam_reader,
		&files.bam_writer);
  bool thread_status = thread.Start();
  if (!thread_status)
    cerr << "threads fail" << endl;
  //StartThreadOrDie(parameters.processors, files.bam_reader);
  

  // free memory and close files
  Deconstruct(&files, &vars);

  return 0;

}



// ====================
// Aux functions
// ====================
void Deconstruct(MainFiles* files, MainVars* vars) {
  // close files
  SR_BamInStreamFree(files->bam_reader);
  bam_close(files->bam_writer);

  // free variables
  SR_BamHeaderFree(vars->bam_header);
}

void InitFiles(const Parameters& parameters, MainFiles* files) {
  // Set the stream mode to "UO" (unique orphan).
  // If the region is specified, than index of the bam is necessary.
  SR_StreamMode streamMode;
  if (parameters.region.empty())
    SR_SetStreamMode(&streamMode, SR_CommonFilter, NULL, SR_NO_SPECIAL_CONTROL);
  else
    // Needs the index of the bam
    SR_SetStreamMode(&streamMode, SR_CommonFilter, NULL, SR_USE_BAM_INDEX);

    
  // Initialize bam input reader.
  // The program will be terminated with printing error message
  // if the given bam cannot be opened.
  files->bam_reader = SR_BamInStreamAlloc(
      parameters.input_bam.c_str(), 
      parameters.mate_window_size,
      parameters.processors,  // number of processors
      2, // the number of alignments can be stored in each chunk of the memory pool
      2, // number of alignments should be cached before report
      &streamMode);

  // Initialize bam output writer
  files->bam_writer = bam_open( parameters.output_bam.c_str(), "w" );
  //files->bam_writer.Open(parameters.output_bam);

  // Initialize reference input reader
  files->ref_reader.open(parameters.input_reference_fasta);
}

void IsInputBamSortedOrDie(const Parameters& parameters,
                           const SR_BamHeader& bam_header) {
  if (!parameters.is_input_sorted &&
      !BamUtilities::IsFileSorted(bam_header.pOrigHeader)) {
    // The input bam is unsorted, exit
    cerr << "ERROR: The input bam seems unsorted. "
         << "Please use bamtools sort to sort the bam" << endl
	 << "       or type -s to ignore this checker." << endl;
	 exit(1);
  }
}

void InitVariablesOrDie(const Parameters& parameters, 
                        MainFiles* files,
                        MainVars* vars) {
  // Load bam header	
  vars->bam_header = SR_BamHeaderAlloc();
  vars->bam_header = SR_BamInStreamLoadHeader(files->bam_reader);

  // Jump the bam if region is given
  if (!parameters.region.empty()) {
    int tid = 0, begin = 0, end = 0;
    if (bam_parse_region(vars->bam_header->pOrigHeader, parameters.region.c_str(), &tid, &begin, &end) == -1) {
      cerr << "ERROR: Parsing the region, " << parameters.region <<", fails." << endl;
      exit(1); // parsing the specified region fails.
    }
    SR_BamInStreamJump(files->bam_reader, tid, begin, end);
  }

  IsInputBamSortedOrDie(parameters, *(vars->bam_header));

  SetAlignmentFilter(parameters, &(vars->alignment_filter));
  SetTargetEvent(parameters, &(vars->target_event));
  SetTargetRegion(parameters, &(vars->target_region));

  // print bam header
  /*
  int n_target = vars->bam_header->pOrigHeader->n_targets;
  cout << vars->bam_header->pOrigHeader->n_targets << endl;
  for (int i = 0; i < n_target; ++i)
    cout << vars->bam_header->pOrigHeader->target_name[i] << endl;
  
  for (int i = 0; i < n_target; ++i)
    cout << vars->bam_header->pOrigHeader->target_len[i] << endl;

  cout << vars->bam_header->pOrigHeader->l_text << endl;
  cout << vars->bam_header->pOrigHeader->n_text << endl;
  cout << vars->bam_header->pOrigHeader->text << endl;

  cout << ((vars->bam_header->pOrigHeader->dict == NULL) ? "NULL" : "Non-null") << endl;
  cout << ((vars->bam_header->pOrigHeader->rg2lib == NULL) ? "NULL" : "Non-null") << endl;
  
  exit(1);
  */
}

void CheckFileOrDie(const Parameters& parameters,
                    const MainFiles& files){

  bool error_found = false;

  if (files.bam_writer == NULL) {
    cerr << "ERROR: Cannot open " 
         << parameters.output_bam 
         << " for writing." << endl
         << "       Please check -o option." << endl;
	  error_found = true;
	}
	
	if (error_found) exit(1);

}

// TODO@WP: Need to check correctness of MD5
void AppendReferenceSequence(bam_header_t* const bam_header, const string& reference_filename){
  FastaReference sp_reader;
  sp_reader.open(reference_filename);

  const int n_special = sp_reader.index->sequenceNames.size();
  char** names   = new char* [n_special];
  uint32_t* lens = new uint32_t [n_special];
  char** md5s    = new char* [n_special];
  const int md5_length = 16;
  for (int i = 0; i < n_special; ++i) {
    string name = sp_reader.index->sequenceNames[i];
    // get and set ref name
    names[i] = new char[name.size() + 1];
    memcpy(names[i], name.c_str(), sizeof(char) * name.size());
    names[i][name.size()] = '\0';
    // get and set ref length
    lens[i] = sp_reader.sequenceLength(name);
    // get ref MD5
    unsigned char md5[md5_length];
    memset(md5, 0, md5_length);
    MD5_CTX context;
    MD5Init(&context);
    char* ptr = new char [lens[i]];
    memcpy(ptr, sp_reader.getSequence(name).c_str(), sizeof(char) * lens[i]);
    MD5Update(&context, (unsigned char*) ptr, lens[i]);
    delete [] ptr;
    MD5Final(md5, &context);
    // set md5
    md5s[i] = new char[32];
    md5s[i][31] = '\0';
    ptr = md5s[i];
    for (int j = 0; j < md5_length; ++j) {
      sprintf(ptr, "%02X", md5[j]);
      ptr += 2;
    }

    //fprintf(stderr,"@SQ\tSN:%s\tLN:%d\tM5:%s\n", names[i], lens[i],md5s[i]);
    //fprintf(stderr,"%s\n", sp_reader.index->sequenceNames[i].c_str());
  }
  BamUtilities::AppendReferenceSequence((const char**)names, lens, (const char**)md5s, n_special, bam_header);

  // delete
  for (int i = 0; i < n_special; ++i) {
    delete [] names[i];
    delete [] md5s[i];
  }
  delete [] names;
  delete [] lens;
  delete [] md5s;
  
  /*
  // open reference file
  FILE* ref_reader = fopen(reference_filename.c_str(), "rb");
  // load the reference header
  int64_t reference_seal;
  SR_RefHeader* reference_header;
  reference_header = SR_RefHeaderRead(&reference_seal, ref_reader);

  if (reference_header->pSpecialRefInfo == NULL) {
    // no special header is detected
    // dose nothing
  } else {
    int n_special = reference_header->pSpecialRefInfo->numRefs;
    int n_normal  = reference_header->numRefs - reference_header->pSpecialRefInfo->numRefs;
    char** names   = new char* [n_special];
    uint32_t* lens = new uint32_t [n_special];
    char** md5s    = new char* [n_special];
    for (int i = 0; i < n_special; ++i) {
      names[i] = reference_header->names[i + n_normal];
      uint32_t begin = SR_SpecialRefGetBeginPos(reference_header, i);
      uint32_t end = reference_header->pSpecialRefInfo->endPos[i];
      lens[i] = end - begin + 1;
      md5s[i] = new char[33];
      memcpy(md5s[i], SR_RefHeaderGetMD5(reference_header, i + n_normal), 32);
      md5s[i][32] = 0;
    }
    BamUtilities::AppendReferenceSequence((const char**)names, lens, (const char**)md5s, n_special, bam_header);
    // delete
    delete [] names;
    delete [] lens;
    for (int i = 0; i < n_special; ++i)
      delete [] md5s[i];
    delete [] md5s;
  }
  
  // close the file
  fclose(ref_reader);
  SR_RefHeaderFree(reference_header);
  */
}

void ResetSoBamHeader(bam_header_t* const bam_header) {
  // Replace SO:coordinate or SO:queryname by SO:unsorted
  BamUtilities::ReplaceHeaderSoText(bam_header);
}

void SetAlignmentFilter(const Parameters& parameters,
                        AlignmentFilter* filter) {
  filter->aligned_base_rate       = parameters.aligned_base_rate;
  filter->allowed_mismatch_rate   = parameters.allowed_mismatch_rate;
  filter->trimming_match_score    = parameters.trimming_match_score;
  filter->trimming_mismatch_score = 0 - parameters.trimming_mismatch_penalty;
  filter->trimming_gap_score      = 0 - parameters.trimming_gap_penalty;
}

void SetTargetEvent(const Parameters& parameters,
                    TargetEvent* target_event) {
  target_event->special_insertion  = parameters.detect_special;
  target_event->medium_sized_indel = !parameters.not_medium_sized_indel;
}

void SetTargetRegion(const Parameters& parameters,
                     TargetRegion* target_region) {
  target_region->fragment_length       = parameters.fragment_length;
  target_region->local_window_size     = parameters.mate_window_size;
  target_region->discovery_window_size = parameters.discovery_window_size;
}
