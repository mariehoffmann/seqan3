// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \brief Provides the configuration for maximum number of errors in percent to the query length across all error types.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/search/algorithm/configuration/utility.hpp>

namespace seqan3::detail
{
/*!\brief A configuration element for the maximum number of errors in percent to the query length across all error types
          (mismatches, insertions, deletions). This is an upper bound of errors independent from error numbers or rates
          of specific error types.
 * \ingroup search_configuration
 */
struct search_config_max_error_rate
{
    //!\brief The actual value.
    std::tuple<double, double, double, double> value;
};

/*!\brief The max_error_rate adaptor enabling pipe notation.
 * \ingroup search_configuration
 */
struct search_config_max_error_rate_adaptor :
    public configuration_fn_base<search_config_max_error_rate_adaptor>
{

    /*!\brief Adds to the configuration a max_error configuration element.
     * \relates seqan3::search_config_max_error
     * \param[in] cfg  The configuration to be extended.
     * \param[in] total_error The number of maximum errors over all error types.
     * \param[in] substitution_error The number of maximum substitution errors.
     * \param[in] insertion_error The number of maximum insertion errors.
     * \param[in] deletion_error The number of maximum deletion errors.
     * \returns A new configuration containing the max_error configuration element.
     */
    template <typename configuration_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    // TODO: replace int by templates and check whether it can be casted to uint8_t
    // TODO: template argument packing. allow any order and subset of strong types
    constexpr auto invoke(configuration_t && cfg, search_cfg::total<double> const total_error,
                                                  search_cfg::substitution<double> const substitution_error,
                                                  search_cfg::insertion<double> const insertion_error,
                                                  search_cfg::deletion<double> const deletion_error) const
    {
        static_assert(is_valid_search_configuration_v<search_cfg::id::max_error, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(search_cfg::id::max_error));

        double total {total_error},
               substitution {substitution_error},
               insertion {insertion_error},
               deletion {deletion_error};

        if (0 > total        || total > 1 ||
            0 > substitution || substitution > 1 ||
            0 > insertion    || insertion > 1 ||
            0 > deletion     || deletion > 1)
            throw std::invalid_argument("Error rates must be between 0 and 1.");

        search_config_max_error_rate tmp{std::tuple<double, double, double, double>{
                                        total,
                                        substitution,
                                        insertion,
                                        deletion}
                                    };
        return std::forward<configuration_t>(cfg).push_front(std::move(tmp));
    }
};

//!\brief Helper template meta-function associated with detail::search_config_max_error_rate.
//!\ingroup search_configuration
template <>
struct on_search_config<search_cfg::id::max_error_rate>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename std::is_same<t, search_config_max_error_rate>::type;
};

//!\brief Mapping from the detail::search_config_max_error_rate type to it's corresponding seqan3::search_cfg::id.
//!\ingroup search_configuration
template <>
struct search_config_type_to_id<search_config_max_error_rate>
{
    //!\brief The associated seqan3::search_cfg::id.
    static constexpr search_cfg::id value = search_cfg::id::max_error_rate;
};
} // namespace seqan3::detail

namespace seqan3::search_cfg
{
/*!\brief A configuration element for the maximum number of errors in percent to the query length across all error types
          (mismatches, insertions, deletions). This is an upper bound of errors independent from error numbers or rates
          of specific error types.
 * \ingroup search_configuration
 */
inline constexpr detail::search_config_max_error_rate_adaptor max_error_rate;

} // namespace seqan3::search_cfg
