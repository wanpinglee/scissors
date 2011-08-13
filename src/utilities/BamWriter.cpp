// ***************************************************************************
// CBamWriter - exports alignment data into the BAM file format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "BamWriter.h"

// constructor
//CBamWriter::CBamWriter(void)
//{}

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

// destructor
CBamWriter::~CBamWriter(void) {
	if( mBGZF.Stream.is_open ) BgzfClose();
}

// closes the BAM file
void CBamWriter::BgzfClose(void) {

	// flush the BGZF block
	BgzfFlushBlock();

	// add an empty block
	mBGZF.BlockOffset = 0;
	int blockLength = BgzfDeflateBlock();
	fwrite(mBGZF.CompressedBlock, 1, blockLength, mBGZF.Stream);

	// flush and close
	fflush( mBGZF.Stream );
	mBGZF.close();
}

// compresses the current block
int CBamWriter::BgzfDeflateBlock(void) {

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
	int compressedLength = 0;

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
void CBamWriter::BgzfFlushBlock(void) {

	// flush all of the remaining blocks
	while(mBGZF.BlockOffset > 0) {

		// compress the data block
		int blockLength = BgzfDeflateBlock();

		// flush the data to our output stream
		int numBytesWritten = fwrite(mBGZF.CompressedBlock, 1, blockLength, mBGZF.Stream);

		if(numBytesWritten != blockLength) {
			printf("ERROR: Expected to write %u bytes during flushing, but wrote %u bytes.\n", blockLength, numBytesWritten);
			exit(1);
		}

		mBGZF.BlockAddress += blockLength;
	}
}

// opens the BAM file for writing
void CBamWriter::BgzfOpen( const string& filename ) {

	mBGZF.Stream.open( filename.c_str(), ofstream::binary );

	if( !mBGZF.Stream.good() ) {
		cout << "ERROR: Unable to open the BAM file " << filename << " for writing." << endl;
		exit( 1 );
	}
}

// writes the supplied data into the BGZF buffer
unsigned int CBamWriter::BgzfWrite(const char* data, const unsigned int dataLen) {

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
void CBamWriter::Close(void) {
	if( mBGZF.Stream.is_open ) BgzfClose();
}


// opens the alignment archive
void CBamWriter::Open(const string& filename, const BamHeader& header) {

	// open the BGZF file for writing
	BgzfOpen( filename );

	// ====================
	// write the SAM header
	// ====================

	// build header tag
	ostringstream sb;
	sb << "@HD\tVN:" << header.Version << "\tSO:";

	switch(header.SortOrder) {
		case SORTORDER_POSITION:
			sb << "coordinate" << endl;
			break;
		case SORTORDER_READNAME:
			sb << "queryname" << endl;
			break;
		default:
			sb << "unsorted" << endl;
	}

	// build the sequence dictionary
	const unsigned int numReferenceSequences = header.pReferenceSequences->size();
	vector<ReferenceSequence>::const_iterator rsIter;
	for(rsIter = header.pReferenceSequences->begin(); rsIter != header.pReferenceSequences->end(); ++rsIter) {
		sb << "@SQ\tSN:" << rsIter->Name << "\tLN:" << rsIter->NumBases;
		if(!rsIter->GenomeAssemblyID.empty()) sb << "\tAS:" << rsIter->GenomeAssemblyID;
		if(!rsIter->MD5.empty())              sb << "\tM5:" << rsIter->MD5;
		if(!rsIter->URI.empty())              sb << "\tUR:" << rsIter->URI;
		if(!rsIter->Species.empty())          sb << "\tSP:" << rsIter->Species;
		sb << endl;
	}

	// build the read groups
	vector<MosaikReadFormat::ReadGroup>::const_iterator rgIter;
	for(rgIter = header.pReadGroups->begin(); rgIter != header.pReadGroups->end(); ++rgIter) {
		sb << "@RG\tID:" << rgIter->ReadGroupID << "\tSM:" << rgIter->SampleName;
		if(!rgIter->LibraryName.empty())      sb << "\tLB:" << rgIter->LibraryName;
		if(!rgIter->Description.empty())      sb << "\tDS:" << rgIter->Description;
		if(!rgIter->PlatformUnit.empty())     sb << "\tPU:" << rgIter->PlatformUnit;
		if(rgIter->MedianFragmentLength != 0) sb << "\tPI:" << rgIter->MedianFragmentLength;
		if(!rgIter->CenterName.empty())       sb << "\tCN:" << rgIter->CenterName;

		switch(rgIter->SequencingTechnology) {
			case ST_454:
				sb << "\tPL:454" << endl;
				break;
			case ST_HELICOS:
				sb << "\tPL:helicos" << endl;
				break;
			case ST_ILLUMINA:
				sb << "\tPL:illumina" << endl;
				break;
			case ST_ILLUMINA_LONG:
				sb << "\tPL:illumina long" << endl;
				break;
			case ST_PACIFIC_BIOSCIENCES:
				sb << "\tPL:pacific biosciences" << endl;
				break;
			case ST_SOLID:
				sb << "\tPL:solid" << endl;
				break;
			case ST_SANGER:
				sb << "\tPL:sanger" << endl;
				break;
			default:
				sb << "\tPL:unknown" << endl;
		}
	}

	if ( !header.pg.ID.empty() ) {
		sb << "@PG\tID:" << header.pg.ID;
		if ( !header.pg.VN.empty() )
			sb << "\tVN:" << header.pg.VN;
		if ( !header.pg.CL.empty() )
			sb << "\tCL:" << header.pg.CL;
		sb << endl;
	}

	// distill the header text
	string samHeader = sb.str();
	sb.str("");

	// ================
	// write the header
	// ================

	// write the BAM signature
	const unsigned char SIGNATURE_LENGTH = 4;
	const char* BAM_SIGNATURE = "BAM\1";
	BgzfWrite(BAM_SIGNATURE, SIGNATURE_LENGTH);

	// write the SAM header text length
	const unsigned int samHeaderLen = samHeader.size();
	BgzfWrite((char*)&samHeaderLen, SIZEOF_INT);

	//printf("samHeaderLen: %u\n%s\n", samHeaderLen, samHeader.c_str());

	// write the SAM header text
	if(samHeaderLen > 0) BgzfWrite(samHeader.data(), samHeaderLen);

	// write the number of reference sequences
	BgzfWrite((char*)&numReferenceSequences, SIZEOF_INT);

	// =============================
	// write the sequence dictionary
	// =============================

	for(rsIter = header.pReferenceSequences->begin(); rsIter != header.pReferenceSequences->end(); ++rsIter) {

		// write the reference sequence name length
		const unsigned int referenceSequenceNameLen = rsIter->Name.size() + 1;
		BgzfWrite((char*)&referenceSequenceNameLen, SIZEOF_INT);

		// write the reference sequence name
		BgzfWrite(rsIter->Name.c_str(), referenceSequenceNameLen);

		// write the reference sequence length
		BgzfWrite((char*)&rsIter->NumBases, SIZEOF_INT);
	}
}

// saves the alignment to the alignment archive
void CBamWriter::SaveAlignment( const bamAlignment& al ) {

	// retrieve our bin
	unsigned int bin = CalculateMinimumBin( al.referenceBegin, al.referenceEnd );

	// assign the BAM core data
	unsigned int buffer[8] = {0};
	buffer[0] = al.referenceIndex;
	buffer[1] = al.referenceBegin;
	buffer[2] = ( bin << 16 ) | ( al.mappingQuality << 8 ) | al.queryName.size();
	buffer[3] = ( al.flag << 16) | al.bamPackedCigar.size();
	buffer[4] = al.sequence.size();

	// mate info
	bool hasPosition     = ( al.referenceBegin != -1 ); 
	bool hasMatePosition = ( al.matePosition != -1 );
	buffer[5] = al.mateReferenceIndex;
	buffer[6] = al.matePosition;
	if ( hasPosition && hasMatePosition 
	&& ( al.referenceIndex == al.mateReferenceIndex ) )
		buffer[7] = al.referenceBegin - al.matePosition;
	else
		buffer[7] = 0;


	// write the block size
	const unsigned int dataBlockSize = al.queryName.size() + al.bamPackedCigar.size() + encodedQueryLen + al.sequence.size();
	const unsigned int blockSize = BAM_CORE_SIZE + dataBlockSize;
	BgzfWrite( (char*) &blockSize, sizeof( int32_t ) );

	// write the BAM core
	BgzfWrite( (char*) &buffer, BAM_CORE_SIZE );

	// write the query name
	BgzfWrite( al.queryName.c_str(), al.queryName.size() );

	// write the packed cigar
	BgzfWrite( al.bamPackedCigar.c_str(), al.bamPackedCigar.size() );

	// write the encoded query sequence
	BgzfWrite(encodedQuery.data(), encodedQueryLen);

	// write the base qualities
	BgzfWrite( al.qual.c_str(), al.qual.size() );
}
