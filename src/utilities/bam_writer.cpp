// ***************************************************************************
// BamWriter - exports alignment data into the BAM file format.
// ***************************************************************************

#include "iostream"
#include "bam_writer.h"

// zlib constants
const uint32_t GZIP_ID1            = 31;
const uint32_t GZIP_ID2            = 139;
const uint32_t CM_DEFLATE          = 8;
const uint32_t FLG_FEXTRA          = 4;
const uint32_t OS_UNKNOWN          = 255;
const uint32_t BGZF_XLEN           = 6;
const uint32_t BGZF_ID1            = 66;
const uint32_t BGZF_ID2            = 67;
const uint32_t BGZF_LEN            = 2;
const uint32_t Z_DEFAULT_MEM_LEVEL = 8;
const int32_t  GZIP_WINDOW_BITS    = -15;


// BGZF block info
const uint32_t BLOCK_HEADER_LENGTH = 18;
const uint32_t BLOCK_FOOTER_LENGTH = 8;

using std::cout;
using std::endl;
using std::min;

// constructor
BamWriter::BamWriter(void) {
}

BamWriter::BamWriter(const string& filename)
	: filename_(filename)
{}

// destructor
BamWriter::~BamWriter(void) {

	if ( outputStream.is_open() )
		Close();
}

// packs an unsigned integer into the specified buffer
inline void BamWriter::BgzfPackUnsignedInt(char* buffer, unsigned int value) {
	buffer[0] = (char)value;
	buffer[1] = (char)(value >> 8);
	buffer[2] = (char)(value >> 16);
	buffer[3] = (char)(value >> 24);
}

// packs an unsigned short into the specified buffer
inline void BamWriter::BgzfPackUnsignedShort(char* buffer, unsigned short value) {
	buffer[0] = (char)value;
	buffer[1] = (char)(value >> 8);
}

// calculates the minimum bin that contains a region [begin, end)
inline unsigned int BamWriter::CalculateMinimumBin(unsigned int begin, unsigned int end) {
	--end;
	if((begin >> 14) == (end >> 14)) return 4681 + (begin >> 14);
	if((begin >> 17) == (end >> 17)) return  585 + (begin >> 17);
	if((begin >> 20) == (end >> 20)) return   73 + (begin >> 20);
	if((begin >> 23) == (end >> 23)) return    9 + (begin >> 23);
	if((begin >> 26) == (end >> 26)) return    1 + (begin >> 26);
	return 0;
}

// compresses the current block
int BamWriter::BgzfDeflateBlock(void) {

	// initialize the gzip header
	char* buffer = mBGZF.CompressedBlock;
	unsigned int bufferSize = mBGZF.CompressedBlockSize;

	memset(buffer, 0, 18);
	buffer[0]  = GZIP_ID1;
	buffer[1]  = (char)GZIP_ID2;
	buffer[2]  = CM_DEFLATE;
	buffer[3]  = FLG_FEXTRA;
	buffer[9]  = (char)OS_UNKNOWN;
	buffer[10] = BGZF_XLEN;
	buffer[12] = BGZF_ID1;
	buffer[13] = BGZF_ID2;
	buffer[14] = BGZF_LEN;

	// loop to retry for blocks that do not compress enough
	int inputLength = mBGZF.BlockOffset;
	unsigned int compressedLength = 0;

	while(true) {

		z_stream zs;
		zs.zalloc    = NULL;
		zs.zfree     = NULL;
		zs.next_in   = (Bytef*)mBGZF.UncompressedBlock;
		zs.avail_in  = inputLength;
		zs.next_out  = (Bytef*)&buffer[BLOCK_HEADER_LENGTH];
		zs.avail_out = bufferSize - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

		// initialize the zlib compression algorithm
		if(deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK) {
			printf("ERROR: zlib deflate initialization failed.\n");
			exit(1);
		}

		// compress the data
		int status = deflate(&zs, Z_FINISH);
		if(status != Z_STREAM_END) {
			deflateEnd(&zs);

			// reduce the input length and try again
			if(status == Z_OK) {
				inputLength -= 1024;
				if(inputLength < 0) {
					printf("ERROR: input reduction failed.\n");
					exit(1);
				}
				continue;
			}

			printf("ERROR: zlib deflate failed.\n");
			exit(1);
		}

		// finalize the compression routine
		if(deflateEnd(&zs) != Z_OK) {
			printf("ERROR: deflate end failed.\n");
			exit(1);
		}

		compressedLength = zs.total_out;
		compressedLength += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

		if(compressedLength > MAX_BGZF_BLOCK_SIZE) {
			printf("ERROR: deflate overflow.\n");
			exit(1);
		}

		break;
	}

	// store the compressed length
	BgzfPackUnsignedShort(&buffer[16], (unsigned short)(compressedLength - 1));

	// store the CRC32 checksum
	unsigned int crc = crc32(0, NULL, 0);
	crc = crc32(crc, (Bytef*)mBGZF.UncompressedBlock, inputLength);
	BgzfPackUnsignedInt(&buffer[compressedLength - 8], crc);
	BgzfPackUnsignedInt(&buffer[compressedLength - 4], inputLength);

	// ensure that we have less than a block of data left
	int remaining = mBGZF.BlockOffset - inputLength;
	if(remaining > 0) {
		if(remaining > inputLength) {
			printf("ERROR: remainder too large.\n");
			exit(1);
		}

		memcpy(mBGZF.UncompressedBlock, mBGZF.UncompressedBlock + inputLength, remaining);
	}

	mBGZF.BlockOffset = remaining;
	return compressedLength;
}

// flushes the data in the BGZF block
void BamWriter::BgzfFlushBlock(void) {

	// flush all of the remaining blocks
	while(mBGZF.BlockOffset > 0) {
		// compress the data block
		int blockLength = BgzfDeflateBlock();

		// write the data to our output stream
		outputStream.write( mBGZF.CompressedBlock, blockLength );
	}
}

// writes the supplied data into the BGZF buffer
unsigned int BamWriter::BgzfWrite(const char* data, const unsigned int dataLen) {

	// initialize
	unsigned int numBytesWritten = 0;
	const char* input = data;
	unsigned int blockLength = mBGZF.UncompressedBlockSize;

	// copy the data to the buffer
	while(numBytesWritten < dataLen) {
		unsigned int copyLength = min(blockLength - mBGZF.BlockOffset, dataLen - numBytesWritten);
		char* buffer = mBGZF.UncompressedBlock;
		memcpy(buffer + mBGZF.BlockOffset, input, copyLength);

		mBGZF.BlockOffset += copyLength;
		input             += copyLength;
		numBytesWritten   += copyLength;

		if(mBGZF.BlockOffset == blockLength) BgzfFlushBlock();
	}

	return numBytesWritten;
}

// closes the alignment archive
void BamWriter::Close(void) {

	// flush the BGZF block
	BgzfFlushBlock();

	// add an empty block
	mBGZF.BlockOffset = 0;
	int blockLength = BgzfDeflateBlock();
	outputStream.write( mBGZF.CompressedBlock, blockLength );
	//fwrite(mBGZF.CompressedBlock, 1, blockLength, mBGZF.Stream );

	// flush and close
	outputStream.flush();
	outputStream.close();
}


// opens the alignment archive
void BamWriter::Open(void) {

	// open the BGZF file for writing
	outputStream.open( filename_.c_str(), ofstream::binary );
	if ( !outputStream.good() ) {
		cout << "ERROR: Unable to open the BAM file " << filename_ << " for writing." << endl;
		exit( 1 );
	}

	// ================
	// write the header
	// ================

	// write the BAM signature
	BgzfWrite( "BAM\1", 4 );

	// write the SAM header text length
	const unsigned int samHeaderLen = 0;
	BgzfWrite( (char*)&samHeaderLen, sizeof( int32_t ) );


	// write the SAM header text
	//if( samHeaderLen > 0 ) BgzfWrite( samHeader.data(), samHeaderLen );

	// write the number of reference sequences
	//BgzfWrite((char*)&numReferenceSequences, SIZEOF_INT);

	// =============================
	// write the sequence dictionary
	// =============================

/*
	for(rsIter = header.pReferenceSequences->begin(); rsIter != header.pReferenceSequences->end(); ++rsIter) {

		// write the reference sequence name length
		const unsigned int referenceSequenceNameLen = rsIter->Name.size() + 1;
		BgzfWrite((char*)&referenceSequenceNameLen, SIZEOF_INT);

		// write the reference sequence name
		BgzfWrite(rsIter->Name.c_str(), referenceSequenceNameLen);

		// write the reference sequence length
		BgzfWrite((char*)&rsIter->NumBases, SIZEOF_INT);
	}
*/

}

// saves the alignment to the alignment archive
void BamWriter::SaveAlignment( const BamAlignment& al ) {

	namespace Constant = BamAlignmentConstant;

	// retrieve our bin
	unsigned int bin = CalculateMinimumBin( al.reference_begin, al.reference_end );

	// assign the BAM core data
	unsigned int buffer[8] = {0};
	buffer[0] = al.reference_index;
	buffer[1] = al.reference_begin;
	buffer[2] = ( bin << 16 ) | ( al.mapping_quality << 8 ) | al.query_name.size();
	buffer[3] = ( al.flag << 16) | al.bam_packed_cigar.size();
	buffer[4] = al.read_length;

	// mate info
	buffer[5] = al.mate_reference_index;
	buffer[6] = al.mate_position;
	buffer[7] = al.isize;

	// write the block size
	const unsigned int dataBlockSize = al.query_name.size() + ( al.bam_packed_cigar.size() * 4 ) + al.encoded_sequence.size() + al.read_length;
	const unsigned int blockSize = Constant::kBamCoreSize + dataBlockSize;
	BgzfWrite( (char*) &blockSize, sizeof( int32_t ) );

	// write the BAM core
	BgzfWrite( (char*) &buffer, Constant::kBamCoreSize );

	// write the query name
	BgzfWrite( al.query_name.c_str(), al.query_name.size() );

	// write the packed cigar
	uint8_t cigar[ al.bam_packed_cigar.size() * 4 ];
	for ( uint32_t i = 0; i < al.bam_packed_cigar.size(); ++i ) {
		uint32_t current = i * 4;
		cigar[ current + 0 ] = al.bam_packed_cigar[i];
		cigar[ current + 1 ] = ( al.bam_packed_cigar[i] >> 8 );
		cigar[ current + 2 ] = ( al.bam_packed_cigar[i] >> 16 );
		cigar[ current + 3 ] = ( al.bam_packed_cigar[i] >> 24 );
	}
	BgzfWrite( (char*)cigar, al.bam_packed_cigar.size() * 4 );

	// write the encoded query sequence
	BgzfWrite( al.encoded_sequence.c_str(), al.encoded_sequence.size() );

	// write the base qualities
	BgzfWrite( al.qual.c_str(), al.qual.size() );
}

