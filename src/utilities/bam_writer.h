// ***************************************************************************
// BamWriter - exports alignment data into the BAM file format.
// ***************************************************************************

#ifndef SRC_UTILITIES_BAMWRITER_H_
#define SRC_UTILITIES_BAMWRITER_H_

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <zlib.h>
#include <stdint.h>
#include <memory>

#include "dataStructures/bam_alignment.h"
#include "samtools/bam.h"

using std::string;
using std::ofstream;
using std::bad_alloc;

// define our sorting types
typedef unsigned char SortOrderType;
const SortOrderType SORTORDER_UNSORTED = 0;
const SortOrderType SORTORDER_READNAME = 10;
const SortOrderType SORTORDER_POSITION = 20;

const uint32_t MAX_BGZF_BLOCK_SIZE = 65536;

// program group in header
struct ProgramGroup {
	string ID;  // program name
	string VN;  // program version
	string CL;  // command line
};

// define our BAM header structure
struct BamHeader {
	SortOrderType SortOrder;
	string Version;
	//vector<> ReadGroup;
	//vector<> pReferenceSequences;
	ProgramGroup pg;

	// constructor
	BamHeader(void)
		: SortOrder(SORTORDER_UNSORTED)
		, Version("1.0")
		//, pReadGroups(NULL)
		//, pReferenceSequences(NULL)
	{}
};

// define our BZGF structure
struct myBGZF {
	unsigned int UncompressedBlockSize;
	unsigned int CompressedBlockSize;
	//unsigned int BlockLength;
	unsigned int BlockOffset;
	uint64_t     BlockAddress;
	//ofstream     Stream;
	char*        UncompressedBlock;
	char*        CompressedBlock;

	// constructor
	myBGZF(void)
		: UncompressedBlockSize(MAX_BGZF_BLOCK_SIZE)
		, CompressedBlockSize(MAX_BGZF_BLOCK_SIZE)
		, BlockOffset(0)
		, BlockAddress(0)
		, UncompressedBlock(NULL)
		, CompressedBlock(NULL)
	{
		try {
			CompressedBlock   = new char[CompressedBlockSize];
			UncompressedBlock = new char[UncompressedBlockSize];
		} catch(bad_alloc) {
			printf("ERROR: Unable to allocate memory for our BGZF object.\n");
			exit(1);
		}
	}

	// destructor
	~myBGZF(void) {
		if(CompressedBlock)   delete [] CompressedBlock;
		if(UncompressedBlock) delete [] UncompressedBlock;
	}

	// copy constructor
	myBGZF ( const myBGZF & copy ) {
		CompressedBlockSize   = copy.CompressedBlockSize;
		UncompressedBlockSize = copy.UncompressedBlockSize;
		CompressedBlock       = new char[ CompressedBlockSize ];
		UncompressedBlock     = new char[ UncompressedBlockSize ];
		memcpy( CompressedBlock, copy.CompressedBlock, CompressedBlockSize );
		memcpy( UncompressedBlock, copy.UncompressedBlock, UncompressedBlockSize );

	}

	// assign operator
	myBGZF& operator=( const myBGZF & copy ) {
		CompressedBlockSize          = copy.CompressedBlockSize;
		UncompressedBlockSize        = copy.UncompressedBlockSize;
		char* temp_CompressedBlock   = new char[ CompressedBlockSize ];
		char* temp_UncompressedBlock = new char[ UncompressedBlockSize ];
		memcpy( temp_CompressedBlock, copy.CompressedBlock, CompressedBlockSize );
		memcpy( temp_UncompressedBlock, copy.UncompressedBlock, UncompressedBlockSize );
		delete [] CompressedBlock;
		delete [] UncompressedBlock;
		CompressedBlock   = temp_CompressedBlock;
		UncompressedBlock = temp_UncompressedBlock;

		return *this;

	}

};

class BamWriter {
public:
	// constructor
	BamWriter(void);
	BamWriter(const string& filename);
	// destructor
	~BamWriter(void);
	// closes the alignment archive
	void Close(void);
	// opens the alignment archive
	void Open(void);
	// saves the alignment to the alignment archive
	void SaveAlignment(const BamAlignment& al);
	void SaveAlignment(const bam1_t& al);
private:
	// compresses the current block
	int BgzfDeflateBlock(void);
	// flushes the data in the BGZF block
	void BgzfFlushBlock(void);
	// packs an unsigned integer into the specified buffer
	static inline void BgzfPackUnsignedInt(char* buffer, unsigned int value);
	// packs an unsigned short into the specified buffer
	static inline void BgzfPackUnsignedShort(char* buffer, unsigned short value);
	// writes the supplied data into the BGZF buffer
	unsigned int BgzfWrite(const char* data, const unsigned int dataLen);
	// calculates the minimum bin that contains a region [begin, end)
	static inline unsigned int CalculateMinimumBin(unsigned int begin, unsigned int end);
	// our BGZF output object
	myBGZF mBGZF;
	ofstream outputStream;
	string filename_;
};

#endif
