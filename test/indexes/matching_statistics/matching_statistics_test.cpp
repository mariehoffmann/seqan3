// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/to_char.hpp>


#include <seqan3/indexes/matching_statistics/all.hpp>
//#include <seqan3/indexes/matching_statistics/unidirectional_matching_statistics.hpp>

#include <gtest/gtest.h>

// forward declaration
//template<typename T> LinkedList<T>::LinkedList()
using namespace seqan3;
using namespace seqan3::literal;

using container_t = typename std::vector<dna4>;
//template <typename container_type> struct seqan3::MS;
namespace fs = std::experimental::filesystem;

//class random_access_iterator_test_fixture : public ::testing::Test
// TODO: inherit from protected seqan3::MS<container_type>
class matching_statistics_test_fixture : public ::testing::Test
{
protected:
    // TODO: use container_t, like dna4_vector instead of std::string
    //container_t s, t;
    std::string s, t;

    virtual void SetUp()
    {
        s = "AACG";     //{dna4::A, dna4::A, dna4::C, dna4::G};
        t = "AACT";     //{dna4::A, dna4::A, dna4::C, dna4::T};
        // code here will execute just before the test ensues
    }
};

TEST(matching_statistics_test, test1)
{
    dna4_vector vec{"ACTTTGATA"_dna4};
    std::string v = vec | view::to_char;
    std::cout << v << std::endl;
}

// unidirectional matching statistics: default constructor
TEST_F(matching_statistics_test_fixture, default_construction)
{
    // default constructors
    MS<container_t> ms{};
    // copy constructor
    MS<container_t> ms2{ms};
    // assignment construction
    MS<container_t> ms3 = ms2;
}

// unidirectional matching statistics: non-default constructors
TEST_F(matching_statistics_test_fixture, non_default_construction)
{
    // construct by sequence
    MS<std::string> ms4(s, t);
    auto paths = ms4.get_output_paths();
    EXPECT_TRUE(fs::exists(paths.first) && fs::exists(paths.second));
    std::string line;
    std::ifstream file1(paths.first);
    if (file1.is_open()) {getline(file1,line); file1.close();}
    EXPECT_EQ(s, line);
    std::ifstream file2(paths.second);
    if (file2.is_open()) {getline(file2,line); file2.close();}
    EXPECT_EQ(t, line);

    // construct by file path
    std::cout << paths.first << std::endl;
    MS<container_t> ms5(paths.first.filename(), paths.second.filename());
    auto paths_ms5 = ms5.get_output_paths();
    std::cout << "ms5: " << paths_ms5.first << std::endl;
    std::cout << fs::exists(paths_ms5.first) << std::endl;
    EXPECT_TRUE(fs::exists(paths_ms5.first) && fs::exists(paths_ms5.second));

}
