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
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sys/types.h>
#include <type_traits>
#include <utility>

#include "range/v3/utility/iterator_traits.hpp"

#include <sdsl/config.hpp> // for cache_config
#include <sdsl/construct.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/util.hpp>

// Note: in seqan3/sdsl-lite/external cmake CMakeLists.txt has to be executed once to create e.g. needed header file divsufsort.h
namespace fs = std::experimental::filesystem;
namespace seqan3
{

// TODO: requires that container's value type fits into a byte! like uint8_t in int_types.hpp
template <typename container_t>
struct MS
{
private:
    // tmp file directory
    fs::path tmp_dir{"./tmp"};
    fs::path full_path_file1{};
    fs::path full_path_file2{};

protected:
    using index_t = signed int; //ranges::v3::size_type_t<container_type>;
    // Algorithm specific value type
    using value_t = unsigned int;
    // Type of container values
    using data_t = ranges::v3::value_type_t<container_t>;
    //!\brief Sequence for which matching statistics will be computed.
    container_t s{NULL};
    container_t t{NULL};

    //!\brief Threshold for number of substring occurrences.
    unsigned short int tau{1};

    sdsl::bit_vector ms;
    sdsl::select_support_mcl<1,1> ss = sdsl::select_support_mcl<1,1>(&ms);

    //!\brief Type of compressed suffix tree.
    typedef sdsl::cst_sada<> cst_t;
    cst_t cst;

    // stream input strings to files
    void write_files()
    {
        assert(s.size() > 0 && t.size() > 0);
        // create directory if empty
        if (!fs::exists(tmp_dir))
            fs::create_directory(tmp_dir);
        // create temporary files if not existent
        if (full_path_file1.empty() || !full_path_file1.has_extension())
            full_path_file1 = tmp_dir / "s.txt";
        if (full_path_file2.empty() || !full_path_file2.has_extension())
            full_path_file2 = tmp_dir / "t.txt";
        // write cached strings to tmp files
        std::ofstream(full_path_file1) << std::string(s.begin(), s.end());
        std::ofstream(full_path_file2) << std::string(t.begin(), t.end());
    }

    //!\brief construct compressed suffix tree of string s according to Sadakane
    void construct_cst()
    {
        assert(s.size() > 0);
        sdsl::tMSS filemap = {{sdsl::conf::KEY_TEXT, infile_text1}, {sdsl::conf::KEY_CST, outfile_cst}};
        sdsl::cache_config config{false, tmp_dir, "", filemap};
        //template <class t_index>
        //void construct(t_index& idx, const std::string& file, cache_config& config, uint8_t num_bytes, cst_tag;
        sdsl::construct(/*t_index&*/ idx, /*const std::string&*/ file, /*cache_config&*/ config, /*uint8_t num_bytes*/ 1, sdsl::cst_tag());

        sdsl::cst_sada(config);

    }


    //!\brief construct Burrows-Wheeler Transform of string s
    void construct_bwt()
    {
        assert(s.size() > 0);
        // bwt needs suffix array (SA) to be in cache. Keys:
        // * conf::KEY_TEXT for t_width=8 or conf::KEY_TEXT_INT for t_width=0
        // * conf::KEY_SA
        // config: Reference to cache configuration
        //cache_config (bool f_delete_files=true, std::string f_dir="./", std::string f_id="", tMSS f_file_map=tMSS())
        sdsl::key_text_trait_impl<8>{s};

        //sdsl::conf::KEY_TEXT = s; // for t_width=8 or conf::KEY_TEXT_INT for t_width=0
        //sdsl::conf::KEY_SA = ;
        sdsl::cache_config config{}; ///*f_delete_files*/ true, /*f_dir*/ "./", /*f_id*/ "", f_file_maptMSS()};

//        sdsl::construct_bwt<8>(&config);
    }

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

    // init with strings and write to file
    constexpr MS(container_t & _s, container_t & _t, unsigned short int const _tau=1) : s(_s), t(_t)
    {
        write_files();
    };

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

    std::pair<fs::path, fs::path> get_output_paths()
    {
        return std::make_pair<fs::path, fs::path>(full_path_file1, full_path_file2);

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
          // build suffix array
          construct_sa();
         // construct BWT from string s
    //     construct_bwt();


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
};

}  // namespace seqan3
