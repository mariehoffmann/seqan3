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

#include <sdsl/suffix_trees>

namespace seqan3
{

template <typename container_t>
struct MS
{
private:
    using index_t =  typename signed int; //ranges::v3::size_type_t<container_type>;
    using value_t = typename unsigned int;
    //!\brief Pointer to sequence for which matching statistics will be computed.
    typename std::add_pointer_t<container_t> s{nullptr};
    //!\brief Pointer to sequence to which matching statistics refer to.
    typename std::add_pointer_t<container_t> t{nullptr};
    //!\brief Threshold for number of substring occurrences.
    unsigned short int tau{1};

    sdsl::bit_vector ms;
    sdsl::select_support_mcl<1,1> ss = sdsl::select_support_mcl<1,1>(&ms);

public:
    //!\brief Default constructor.
    constexpr MS() = default;

    //!\brief Copy constructor.
    constexpr MS(MS const &) = default;

    //!\brief Copy construction via assignment.
    constexpr MS & operator=(MS const &) = default;

    //!\brief Move constructor.
    constexpr MS (MS &&) = default;

    //!\brief Move assignment.
    constexpr MS & operator=(MS &&) = default;

    constexpr MS(container_t & _s, container_t & _t) : s(_s), t(_t) {}
    constexpr MS(container_t & _s, container_t & _t, , unsigned short int const _tau) : s(_s), t(_t), tau{_tau} {}

    //!\brief Use default deconstructor.
    ~MS() = default;

    constexpr value_t select(index_t const i)
    {
        assert(i >= 0);
        return ss(i);
    }

    constexpr value_t operator[](index_t const i) const
    {
        assert(i >= -1);
        if (i == -1)
            return 1;
        return select(i) - 2*i;
    }

    /*
     * Definition "unidirectional matching statistics":
     * Given two strings s and t and a threshold tau > 0,the unidirectional matching
     * statistics MS(t,s,tau) of t with respect to s is a vector of length |t| that
     * stores at index i in [0..|t| − 1] the length of the longest prefix of
     * t[i..|t| − 1] that occurs at least tau times in s.
     */
    void compute()
    {

    }
}

template <typename container_type>
struct DS
{

private:
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



}  // namespace seqan3
