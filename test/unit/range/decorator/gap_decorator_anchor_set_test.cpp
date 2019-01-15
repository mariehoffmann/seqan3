// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/range/decorator/gap_decorator_anchor_set.hpp>

#include "./gap_decorator_test_template.hpp"
#include "../../alignment/aligned_sequence_test_template.hpp"

using namespace seqan3;

typedef typename std::vector<dna4> sequence_type;

class aligned_sequence<gap_decorator_anchor_set<sequence_type>> : public ::testing::Test
{
public:
    // Initializer function is needed for the typed test because the gapped_decorator
    // will be initialized differently than the naive vector<gapped<dna>>.
    void initialize_typed_test_container(gap_decorator_anchor_set<sequence_type> & container, std::string const && target)
    {
        for (auto & val : target)
        {
            typename sequence_type::value_type cval{};
            assign_char(cval, val);
            container.push_back(cval);
        }
    }
};

using gap_decorator_types = ::testing::Types<gap_decorator_anchor_set<sequence_type>>;

INSTANTIATE_TYPED_TEST_CASE_P(gap_decorator, gap_decorator, gap_decorator_types);
INSTANTIATE_TYPED_TEST_CASE_P(aligned_sequence, aligned_sequence, gap_decorator_types);