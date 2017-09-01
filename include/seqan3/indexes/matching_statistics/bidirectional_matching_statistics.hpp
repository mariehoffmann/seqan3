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
 * \brief Bidirectional matching statistics (BMS) for sequence comparison.
 * \ingroup container
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <cassert>

#include <sdsl/suffix_trees.hpp>
#include <sdsl/construct_bwt.hpp>

namespace seqan3
{

// TODO: requires that container's value type fits into a byte! like uint8_t in int_types.hpp
template <typename container_t>
struct BMS : MS
{
private:
protected:
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

    void construct_bwt()
    {
        assert(s.size() > 0)
        sdsl::construct_bwt<1>(cache_config& config)
    }

     /*
      * Definition "unidirectional matching statistics":
      * Given two strings s and t and a threshold tau > 0,the unidirectional matching
      * statistics MS(t,s,tau) of t with respect to s is a vector of length |t| that
      * stores at index i in [0..|t| − 1] the length of the longest prefix of
      * t[i..|t| − 1] that occurs at least tau times in s.
      * t = AACT, s = AACG, MS = 3210, ms = 000111
      */
      void compute()
      {
         // construct BWT from string s
         construct_bwt();


        ms.resize(2*t.size());
        index_t tmp = -1; // MS[-1]
        index_t j = 0;
        // 1st pass: compute consecutive 1s for ms
        sdsl::bit_vector runs = sdsl::bit_vector(t.size()-1, 0);
        for (index_t i = 0; i < t.size(); ++i)
        {

            if (j >= ms.size()){
                std::cout << "error, j is exceeding allocated bit_vector size with j = " << j << std::endl;
                break;
            }
        }

    }
}

}  // namespace seqan3
