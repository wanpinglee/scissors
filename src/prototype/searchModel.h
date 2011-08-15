
#ifndef _SEARCHMODEL_H_
#define _SEARCHMODEL_H_


// define the search model
struct SearchModel {
	
	// NOTE: Is there any other model?
	enum targetSequenceModel {
		bool forward;
		bool reverseComplement;
		bool inverse;
	};
	
	bool IsTargetBeforeAnchorMate;
	bool IsTargetReverseComplement;
	targetSequenceModel sequenceModel;

	searchModel( void )
		: IsTargetBeforeAnchorMate(false)
		, IsTargetReverseComplement(false)
		, sequenceModel(forward)
	{}
};

#endif
