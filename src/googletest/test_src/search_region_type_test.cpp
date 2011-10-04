#include "dataStructures/search_region_type.h"

#include "gtest/gtest.h"

namespace {

// for illumina forward anchors
const static SearchRegionType::RegionType kRegionType1 = {true, true, true};
const static SearchRegionType::RegionType kRegionType2 = {false, true, true};
const static SearchRegionType::RegionType kRegionType3 = {true, false, true};
const static SearchRegionType::RegionType kRegionType4 = {false, false, true};

// for illumina reversed complement anchors
const static SearchRegionType::RegionType kRegionType5 = {false, false, false};
const static SearchRegionType::RegionType kRegionType6 = {true, false, false};
const static SearchRegionType::RegionType kRegionType7 = {false, true, false};
const static SearchRegionType::RegionType kRegionType8 = {true, true, false};

void Validate( const RegionType& expect, const RegionType& actual) {
  EXPECT_EQ(expect.upstream, actual.upstream);
  EXPECT_EQ(expect.sequence_inverse, actual.sequence_inverse);
  EXPECT_EQ(expect.sequence_complement, actual.sequence_complement);
}

void TestLoading(const bool strand) {
  SearchRegionType::RegionType region_type;
  SearchRegionType search_region_type;  // the default tech is ILLUMINA
  bool region_type_obtained = false;
  
  search_region_type.ResetRegionTypeList();

  // get first type
  region_type_obtained = search_region_type.GetNextRegionType(strand, &region_type);
  EXPECT_TRUE(region_type_obtained);
  if (strand)
    Validate(kRegionType1, region_type);
  else
    Validate(kRegionType5, region_type);

  // get second type
  region_type_obtained = search_region_type.GetNextRegionType(strand, &region_type);
  EXPECT_TRUE(region_type_obtained);
  if (strand)
    Validate(kRegionType2, region_type);
  else
    Validate(kRegionType6, region_type);

  // get third type
  region_type_obtained = search_region_type.GetNextRegionType(strand, &region_type);
  EXPECT_TRUE(region_type_obtained);
  if (strand)
    Validate(kRegionType3, region_type);
  else
    Validate(kRegionType7, region_type);

  // get fourth type
  region_type_obtained = search_region_type.GetNextRegionType(strand, &region_type);
  EXPECT_TRUE(region_type_obtained);
  if (strand)
    Validate(kRegionType4, region_type);
  else
    Validate(kRegionType8, region_type);

  // try to get more
  region_type_obtained = search_region_type.GetNextRegionType(strand, &region_type);
  EXPECT_FALSE(region_type_obtained);

  EXPECT_FALSE(search_region_type.SetCurrentTypeSuccess(strand));
  EXPECT_FALSE(search_region_type.SetCurrentTypeSuccess(!strand));

  // rewind the list
  search_region_type.RewindRegionTypeList();

  // get first type
  region_type_obtained = search_region_type.GetNextRegionType(!strand, &region_type);
  EXPECT_TRUE(region_type_obtained);
  if (strand)
    Validate(kRegionType5, region_type);
  else
    Validate(kRegionType1, region_type);

  // get second type
  region_type_obtained = search_region_type.GetNextRegionType(!strand, &region_type);
  EXPECT_TRUE(region_type_obtained);
  if (strand)
    Validate(kRegionType6, region_type);
  else
    Validate(kRegionType2, region_type);
}

TEST(SearchRegionTypeTest, Loading){
  TestLoading(true);
  TestLoading(false);
}


void TestSetting(const bool strand) { 
  SearchRegionType::RegionType region_type;
  SearchRegionType search_region_type;  // the default tech is ILLUMINA

  bool is_set = search_region_type.SetCurrentTypeSuccess(strand);
  EXPECT_FALSE(is_set);

  // rewind the list
  search_region_type.RewindRegionTypeList();
  is_set = search_region_type.SetCurrentTypeSuccess(strand);
  EXPECT_FALSE(is_set);

  // get first type
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType1, region_type);
  else
    Validate(kRegionType5, region_type);

  // get second type
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType2, region_type);
  else
    Validate(kRegionType6, region_type);

  is_set = search_region_type.SetCurrentTypeSuccess(strand);
  EXPECT_TRUE(is_set);

  search_region_type.RewindRegionTypeList();

  // get first type
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType2, region_type);
  else
    Validate(kRegionType6, region_type);

  // get second type
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType1, region_type);
  else
    Validate(kRegionType5, region_type);

  // get third type
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType3, region_type);
  else
    Validate(kRegionType7, region_type);

  // get fourth type
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType4, region_type);
  else
    Validate(kRegionType8, region_type);

  is_set = search_region_type.SetCurrentTypeSuccess(strand);
  EXPECT_TRUE(is_set);

  search_region_type.RewindRegionTypeList();

  search_region_type.GetNextRegionType(strand, &region_type);
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType1, region_type);
  else
    Validate(kRegionType5, region_type);
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType4, region_type);
  else
    Validate(kRegionType8, region_type);
  search_region_type.GetNextRegionType(strand, &region_type);
  if (strand)
    Validate(kRegionType3, region_type);
  else
    Validate(kRegionType7, region_type);

}

} // namespace
