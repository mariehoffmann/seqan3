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
//#ifndef __MS_HPP
//#define __MS_HPP

#include <algorithm>
#include <cassert>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/types.h>
#include <type_traits>
#include <utility>
#include <vector>

#include <range/v3/utility/iterator_traits.hpp>
#include <range/v3/view/reverse.hpp>

#include <sdsl/config.hpp> // for cache_config
#include <sdsl/construct.hpp>
#include <sdsl/suffix_trees.hpp>
//#include <sdsl/util.hpp>

//#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
//#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/to_char.hpp>

// Note: in seqan3/sdsl-lite/external cmake CMakeLists.txt has to be executed once to create e.g. needed header file divsufsort.h
namespace fs = std::experimental::filesystem;

using namespace seqan3;
using namespace seqan3::literal;

#define STRING_S 0
#define STRING_T 1

//template <typename alphabet_type>
//auto const seqan3::view::char_to<alphabet_type>;

// TODO: requires that container's value type fits into a byte! like uint8_t in int_types.hpp
template <typename container_t=std::string>
struct MS
{
private:

protected:
    using index_t = signed int; //ranges::v3::size_type_t<container_type>;
    // Algorithm specific value type
    using value_t = unsigned int;
    // Type of container values
    using alphabet_t = ranges::v3::value_type_t<container_t>;

    //!\brief Threshold for number of substring occurrences.
    unsigned short int tau{1};

    sdsl::bit_vector ms;
    sdsl::select_support_mcl<1,1> ss = sdsl::select_support_mcl<1,1>(&ms);

    //!\brief Type of compressed suffix tree.
    typedef sdsl::cst_sada<> cst_t;

    struct source_t
    {
        fs::path filename;
        fs::path filename_rev;
        typename std::add_pointer_t<container_t> sequence{nullptr};
    };

public:
    // tmp file directory
    fs::path tmp_dir{"./tmp"};

    std::array<source_t, 2> srcs = {source_t{fs::path{}}, source_t{fs::path{}}};
    //cst_t cst;  // TODO: move to private section, derive test class from this for access

    // TODO: use alphabet vectors and crate strings from them
    // stream input strings to files

    auto get_absolute_path(fs::path filename)
    {
        return fs::current_path() / tmp_dir / filename;
    }

    auto get_output_paths()
    {
        return std::make_pair(get_absolute_path(srcs[0].filename), get_absolute_path(srcs[1].filename));
    }

    auto get_strings()
    {
        return std::make_pair(*srcs[0]->sequence, *srcs[1]->sequence);
    }

    void write_files()
    {
        assert(srcs[0].sequence->size() > 0 && srcs[1].sequence->size() > 0);
        // create reverse strings
//        srcs[STRING_S_REV].sequence = &std::string(srcs[STRING_S].sequence->rbegin(), srcs[STRING_S].sequence->rend());
        // create directory if empty
        if (!fs::exists(tmp_dir))
            fs::create_directory(tmp_dir);
        // create temporary files if not existent
        if (srcs[0].filename.empty() || !srcs[0].filename.has_extension()){
            srcs[0].filename = "s.txt";
            srcs[0].filename_rev = "s_rev.txt";
        }
        if (srcs[1].filename.empty() || !srcs[1].filename.has_extension())
        {
            srcs[1].filename = "t.txt";
            srcs[1].filename_rev = "t_rev.txt";
        }
        //auto v = ssss | view::char_to<dna4>;
        //std::cout << v << std::endl;
        // write cached strings and their reverse into files
        for (auto i = 0; i < srcs.size(); ++i)
        {
            std::ofstream(get_absolute_path(srcs[i].filename)) << *(srcs[i].sequence);
            std::cout << "get_absolute_path(srcs[i].filename) exists: " << fs::exists(get_absolute_path(srcs[i].filename)) << std::endl;
            std::string tmp(srcs[i].sequence->rbegin(), srcs[i].sequence->rend());
            std::ofstream(get_absolute_path(srcs[i].filename_rev)) << tmp;
            std::cout << "get_absolute_path(srcs[i].filename_rev) exists: " << fs::exists(get_absolute_path(srcs[i].filename_rev)) << std::endl;
        }
    }

    //!\brief construct compressed suffix tree of string s according to Sadakane
    void construct_cst(size_t const i, cst_t & cst, bool const reverse=false)
    {
        fs::path src_file = (reverse ? srcs[i].filename_rev : srcs[i].filename);
        assert(fs::exists(get_absolute_path(srcs[i].filename)));
        sdsl::construct(cst, get_absolute_path(srcs[i].filename), 1);

        auto root = cst.root();
        for (auto child: cst.children(root)) {
            std::cout << "sada id = " << cst.id(child) << std::endl;
            if (cst.id(child) > 0)
                std::cout << "1st char of edge = '" << cst.edge(child, 1) << "'" << std::endl;
        }

        //sdsl::tMSS filemap = {{sdsl::conf::KEY_TEXT, infile_text1}, {sdsl::conf::KEY_CST, outfile_cst}};
        //sdsl::cache_config config{false, tmp_dir, "", filemap};
        //template <class t_index>
        //void construct(t_index& idx, const std::string& file, cache_config& config, uint8_t num_bytes, cst_tag;
        //sdsl::construct(/*t_index&*/ idx, /*const std::string&*/ file, /*cache_config&*/ config, /*uint8_t num_bytes*/ 1, sdsl::cst_tag());

        //sdsl::cst_sada(config);

    }


    //!\brief construct Burrows-Wheeler Transform of string s
    // TODO: bug when using "std::array<source_t, 2>::size_type" instead of size_t
    void construct_bwt(size_t i)
    {
        assert(srcs[i].sequence->size() > 0);
        // bwt needs suffix array (SA) to be in cache. Keys:
        // * conf::KEY_TEXT for t_width=8 or conf::KEY_TEXT_INT for t_width=0
        // * conf::KEY_SA
        // config: Reference to cache configuration
        //cache_config (bool f_delete_files=true, std::string f_dir="./", std::string f_id="", tMSS f_file_map=tMSS())
        //sdsl::key_text_trait_impl<8>{s};

        //sdsl::conf::KEY_TEXT = s; // for t_width=8 or conf::KEY_TEXT_INT for t_width=0
        //sdsl::conf::KEY_SA = ;
        sdsl::cache_config config{}; ///*f_delete_files*/ true, /*f_dir*/ "./", /*f_id*/ "", f_file_maptMSS()};

//        sdsl::construct_bwt<8>(&config);
    }


    //!\brief Default constructor.
    MS() = default;

    //!\brief Copy constructor.
    constexpr MS(MS const &) = default;

    //!\brief Copy construction via assignment.
    constexpr MS & operator=(MS const &) = default;

    //!\brief Move constructor.
    constexpr MS (MS &&) = default;

    //!\brief Move assignment.
    constexpr MS & operator=(MS &&) = default;

    // init with strings and write to file
    //explicit constexpr random_access_iterator(container_type & host) noexcept : host{&host} {}
    explicit constexpr MS(container_t & _s, container_t & _t, unsigned short int const _tau=1) noexcept
    {
        srcs[STRING_S].sequence = &_s;
        srcs[STRING_T].sequence = &_t;
        write_files();
    }

    // _filex has to be located in current_abs_dir/tmp_dir/filename, otherwise assertion is thrown
    constexpr MS(fs::path _file1, fs::path _file2) {
        srcs[0].filename = _file1;
        srcs[1].filename = _file2;
        assert(fs::exists(get_absolute_path(_file1)) && fs::exists(get_absolute_path(_file1)));
    }

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
      * t = AACT, s = AACG, MS = 3210, ms = 000111
      */
      void compute()
      {
         // build suffix tree (ST) from string s
         cst_t cst_s, cst_s_rev;
         construct_cst(0, cst_s);
         // TODO: named attribute not possible like reverse=true
         construct_cst(0, cst_s_rev, true);

         // construct Burrows-Wheeler transform (BWT) from string s
//         construct_bwt(0);


        ms.resize(2*srcs[1].sequence->size());
        index_t tmp = -1; // MS[-1]
        index_t j = 0;
        // 1st pass: compute consecutive 1s for ms
        /*
        sdsl::bit_vector runs = sdsl::bit_vector(srcs[1].sequence->size()-1, 0);
        for (index_t i = 0; i < srcs[1].sequence->size(); ++i)
        {

            if (j >= ms.size()){
                std::cout << "error, j is exceeding allocated bit_vector size with j = " << j << std::endl;
                break;
            }
        }*/
    }
};

//#endif

//}  // namespace seqan3
