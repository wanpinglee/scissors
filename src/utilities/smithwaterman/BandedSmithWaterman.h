#ifndef _BANDEDSMITHWATERMAN_H_
#define _BANDEDSMITHWATERMAN_H_

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <string.h>
#include <stdio.h>

#include "dataStructures/alignment.h"

using namespace std;


struct BandedSmithWatermanHashRegion {
	uint32_t reference_begin;
	//uint32_t reference_end;   // is not used in current design 
	uint32_t query_begin;
	//uint32_t query_end;       // is not used in current design

	BandedSmithWatermanHashRegion()
		: reference_begin(0)
	//	, reference_end(0)
		, query_begin(0)
	//	, query_end(0)
	{}

};

class CBandedSmithWaterman {
public:
	// constructor
	CBandedSmithWaterman();
	CBandedSmithWaterman(float matchScore, float mismatchScore, 
		float gapOpenPenalty, float gapExtendPenalty, unsigned int bandWidth);
	// destructor
	~CBandedSmithWaterman(void);
	// aligns the query sequence to the anchor using the Smith Waterman Gotoh algorithm
	void Align 
		( Alignment& al
		, const char* s1
		, const unsigned int& s1Length
		, const char* s2
		, const unsigned int& s2Length
		, const BandedSmithWatermanHashRegion& hash_region);
	// enables homo-polymer scoring
	void EnableHomoPolymerGapPenalty(float hpGapOpenPenalty);

private:
	CBandedSmithWaterman ( const CBandedSmithWaterman& copy );
	CBandedSmithWaterman& operator= ( const CBandedSmithWaterman& copy );

	// ====
	// data
	// ====
	enum PositionType {
		Position_REF_AND_QUERY_ZERO,
		Position_REF_ZERO,
		Position_QUERY_ZERO,
		Position_REF_AND_QUERO_NONZERO
	};

	enum DirectionType {
		Directions_STOP,
		Directions_LEFT,
		Directions_DIAGONAL,
		Directions_UP
	};

	struct ElementInfo {
		unsigned int Direction             : 2;
		unsigned int mSizeOfVerticalGaps   : 15;
		unsigned int mSizeOfHorizontalGaps : 15;
	};

	// =========
	// functions
	// =========
        void Init();	
	// calculates the score during the forward algorithm
	float CalculateScore
		( const char* s1
		, const char* s2
		, const unsigned int rowNum
		, const unsigned int columnNum
		, float& currentQueryGapScore
		, const unsigned int rowOffset
		, const unsigned int columnOffset);
	// creates a simple scoring matrix to align the nucleotides and the ambiguity code N
	void CreateScoringMatrix ( void );
	// corrects the homopolymer gap order for forward alignments
	void CorrectHomopolymerGapOrder ( const unsigned int numBases, const unsigned int numMismatches );
	// returns the maximum floating point number
	static inline float MaxFloats ( const float& a, const float& b, const float& c );
	// reinitializes the matrices
	void ReinitializeMatrices 
		( const PositionType& positionType
		, const unsigned int& s1Length
		, const unsigned int& s2Length
		, const BandedSmithWatermanHashRegion& hr);
	// performs the backtrace algorithm
	void Traceback 
		( Alignment& alignment
		, const char* s1
		, const char* s2
		, const unsigned int s2Length
		, unsigned int bestRow
		, unsigned int bestColumn
		, const unsigned int rowOffset
		, const unsigned int columnOffset );
	// updates the best score during the forward algorithm
	inline void UpdateBestScore 
		( unsigned int& bestRow
		, unsigned int& bestColumn
		, float& bestScore
		, const unsigned int rowNum
		, const unsigned int columnNum
		, const float score);
	// our simple scoring matrix
	float mScoringMatrix[MOSAIK_NUM_NUCLEOTIDES][MOSAIK_NUM_NUCLEOTIDES];
	// keep track of maximum initialized sizes
	unsigned int mCurrentMatrixSize;
	unsigned int mCurrentAnchorSize;
	unsigned int mCurrentAQSumSize;
	unsigned int mBandwidth;
	// store the backtrace pointers
	ElementInfo* mPointers;
	// define scoring constants
	const float mMatchScore;
	const float mMismatchScore;
	const float mGapOpenPenalty;
	const float mGapExtendPenalty;
	// score if xi aligns to a gap after yi
	float* mAnchorGapScores;
	// best score of alignment x1...xi to y1...yi
	float* mBestScores;
	// our reversed alignment
	char* mReversedAnchor;
	char* mReversedQuery;
	// toggles the use of the homo-polymer gap open penalty
	bool  mUseHomoPolymerGapOpenPenalty;
	float mHomoPolymerGapOpenPenalty;
};

#endif
