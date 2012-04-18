#include "gtest/gtest.h"

#include "dataStructures/alignment.h"
#include "utilities/miscellaneous/alignment_filter.h"

TEST(AlignmentFilter, NoTrim) {
  Alignment al;
  al.reference = "ACTGACTGACGTACGTA";
  al.query     = "ACTGACTGACGTACGTA";
  al.reference_begin = 0;
  al.reference_end   = 16;
  al.query_begin     = 0;
  al.query_end       = 16;
  al.num_mismatches  = 0;

  Alignment golden = al;

  AlignmentFilter filter;
  ASSERT_FALSE(AlignmentFilterApplication::TrimAlignment(filter, &al));
  
  ASSERT_EQ(golden.reference_begin, al.reference_begin);
  ASSERT_EQ(golden.reference_end, al.reference_end);
  ASSERT_EQ(golden.query_begin, al.query_begin);
  ASSERT_EQ(golden.query_end, al.query_end);
  ASSERT_EQ(golden.num_mismatches, al.num_mismatches);
  ASSERT_TRUE(golden.reference == al.reference);
  ASSERT_TRUE(golden.query == al.query);
}

TEST(AlignmentFilter, TrimBegin) {
  Alignment al;
  al.reference = "ACTGACTGACGTACGTA";
  al.query     = "GGTGACTGACGTACGTA";
  al.reference_begin = 0;
  al.reference_end   = 16;
  al.query_begin     = 0;
  al.query_end       = 16;
  al.num_mismatches  = 2;

  Alignment golden;
  golden.reference = "TGACTGACGTACGTA";
  golden.query     = "TGACTGACGTACGTA";
  golden.reference_begin = 2;
  golden.reference_end   = 16;
  golden.query_begin     = 2;
  golden.query_end       = 16;
  golden.num_mismatches  = 0;


  AlignmentFilter filter;
  ASSERT_TRUE(AlignmentFilterApplication::TrimAlignment(filter, &al));
  
  ASSERT_EQ(golden.reference_begin, al.reference_begin);
  ASSERT_EQ(golden.reference_end, al.reference_end);
  ASSERT_EQ(golden.query_begin, al.query_begin);
  ASSERT_EQ(golden.query_end, al.query_end);
  ASSERT_EQ(golden.num_mismatches, al.num_mismatches);
  ASSERT_TRUE(golden.reference == al.reference);
  ASSERT_TRUE(golden.query == al.query);
}

TEST(AlignmentFilter, TrimBeginGap) {
  Alignment al;
  al.reference = "ACTGACTGACGTACGTA";
  al.query     = "A--GACTGACGTACGTA";
  al.reference_begin = 0;
  al.reference_end   = 16;
  al.query_begin     = 0;
  al.query_end       = 14;
  al.num_mismatches  = 2;

  Alignment golden;
  golden.reference = "GACTGACGTACGTA";
  golden.query     = "GACTGACGTACGTA";
  golden.reference_begin = 3;
  golden.reference_end   = 16;
  golden.query_begin     = 1;
  golden.query_end       = 14;
  golden.num_mismatches  = 0;


  AlignmentFilter filter;
  ASSERT_TRUE(AlignmentFilterApplication::TrimAlignment(filter, &al));
  
  ASSERT_EQ(golden.reference_begin, al.reference_begin);
  ASSERT_EQ(golden.reference_end, al.reference_end);
  ASSERT_EQ(golden.query_begin, al.query_begin);
  ASSERT_EQ(golden.query_end, al.query_end);
  ASSERT_EQ(golden.num_mismatches, al.num_mismatches);
  ASSERT_TRUE(golden.reference == al.reference);
  ASSERT_TRUE(golden.query == al.query);
}

TEST(AlignmentFilter, TrimEnd) {
  Alignment al;
  //                             X
  al.reference = "ACTGACTGACGTACGTA";
  al.query     = "ACTGACTGACGTACGGA";
  al.reference_begin = 0;
  al.reference_end   = 16;
  al.query_begin     = 0;
  al.query_end       = 16;
  al.num_mismatches  = 1;

  Alignment golden;
  golden.reference = "ACTGACTGACGTACG";
  golden.query     = "ACTGACTGACGTACG";
  golden.reference_begin = 0;
  golden.reference_end   = 14;
  golden.query_begin     = 0;
  golden.query_end       = 14;
  golden.num_mismatches  = 0;


  AlignmentFilter filter;
  ASSERT_TRUE(AlignmentFilterApplication::TrimAlignment(filter, &al));
  
  ASSERT_EQ(golden.reference_begin, al.reference_begin);
  ASSERT_EQ(golden.reference_end, al.reference_end);
  ASSERT_EQ(golden.query_begin, al.query_begin);
  ASSERT_EQ(golden.query_end, al.query_end);
  ASSERT_EQ(golden.num_mismatches, al.num_mismatches);
  ASSERT_TRUE(golden.reference == al.reference);
  ASSERT_TRUE(golden.query == al.query);
}

TEST(AlignmentFilter, TrimEndGap) {
  Alignment al;
  al.reference = "ACTGACTGACGTACGTA";
  al.query     = "ACTGACTGACGTACG-A";
  al.reference_begin = 0;
  al.reference_end   = 16;
  al.query_begin     = 0;
  al.query_end       = 15;
  al.num_mismatches  = 1;

  Alignment golden;
  golden.reference = "ACTGACTGACGTACG";
  golden.query     = "ACTGACTGACGTACG";
  golden.reference_begin = 0;
  golden.reference_end   = 14;
  golden.query_begin     = 0;
  golden.query_end       = 14;
  golden.num_mismatches  = 0;


  AlignmentFilter filter;
  ASSERT_TRUE(AlignmentFilterApplication::TrimAlignment(filter, &al));
  
  ASSERT_EQ(golden.reference_begin, al.reference_begin);
  ASSERT_EQ(golden.reference_end, al.reference_end);
  ASSERT_EQ(golden.query_begin, al.query_begin);
  ASSERT_EQ(golden.query_end, al.query_end);
  ASSERT_EQ(golden.num_mismatches, al.num_mismatches);
  ASSERT_TRUE(golden.reference == al.reference);
  ASSERT_TRUE(golden.query == al.query);
}

TEST(AlignmentFilter, TrimBoth) {
  Alignment al;
  al.reference = "ACTGACTGACGTACGTA";
  al.query     = "ACCGACTGGCGTACG-A";
  al.reference_begin = 0;
  al.reference_end   = 16;
  al.query_begin     = 0;
  al.query_end       = 15;
  al.num_mismatches  = 2;

  Alignment golden;
  golden.reference = "GACTGACGTACG";
  golden.query     = "GACTGGCGTACG";
  golden.reference_begin = 3;
  golden.reference_end   = 14;
  golden.query_begin     = 3;
  golden.query_end       = 14;
  golden.num_mismatches  = 0;


  AlignmentFilter filter;
  ASSERT_TRUE(AlignmentFilterApplication::TrimAlignment(filter, &al));
  
  ASSERT_EQ(golden.reference_begin, al.reference_begin);
  ASSERT_EQ(golden.reference_end, al.reference_end);
  ASSERT_EQ(golden.query_begin, al.query_begin);
  ASSERT_EQ(golden.query_end, al.query_end);
  ASSERT_EQ(golden.num_mismatches, al.num_mismatches);
  ASSERT_TRUE(golden.reference == al.reference);
  ASSERT_TRUE(golden.query == al.query);
}
