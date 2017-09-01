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

/*!\file
 * \brief Matching statistics (e.g. shortest unique substrings) for sequence comparison.
 * \ingroup container
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <cassert>

#include <sdsl/suffix_trees.hpp>
#include <sdsl/construct_bwt.hpp>

namespace seqan3
{

template <typename container_type>
struct DS
{
protected:
    using position_type =  ranges::v3::size_type_t<container_type>;
    //!\brief Pointer to sequence for which unidirectional distinguishing statistics will be computed.
    typename std::add_pointer_t<container_t> t{nullptr};
    //!\brief Threshold for number of substring occurrences.
    unsigned short int tau{1};

public:
    //!\brief Default constructor.
    constexpr DS() = default;

    //!\brief Copy constructor.
    constexpr DS(DS const &) = default;

    //!\brief Copy construction via assignment.
    constexpr DS & operator=(DS const &) = default;

    //!\brief Move constructor.
    constexpr DS (DS &&) = default;

    //!\brief Move assignment.
    constexpr DS & operator=(DS &&) = default;

    constexpr DS(container_t & _t) : t(_t) {}
    constexpr DS(container_t & _t, unsigned short int const _tau) : t(_t), tau{_tau} {}

    //!\brief Use default deconstructor.
    ~DS() = default;

    /*
     * Definition "unidirectional distinguishing statistics":
     * Given a string t and a threshold tau > 0,the unidirectional distinguishing
     * statistics DS(t,tau) of t is a vector of length |t| that stores at index
     * i in [0..|t|−1] the length of the shortest prefix of t[i..|t|−1]$ that occurs
     * at most tau times in t.
     */
    void compute()
    {

    }

}

template <typename container_type>
struct BMS
{

}

template <typename container_type>
struct BDS
{

}

}
