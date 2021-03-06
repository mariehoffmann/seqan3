// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::all.
 */

#pragma once

#include <string_view>

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/detail.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

namespace seqan3::detail
{

// ============================================================================
//  all_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for views::all.
 */
class all_fn : public adaptor_base<all_fn>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = adaptor_base<all_fn>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief Type erase if possible and delegate to std::views::all otherwise.
     * \returns An instance of std::span, std::basic_string_view, std::ranges::subrange or std::views::all's return type.
     */
    template <std::ranges::range urng_t>
    static constexpr auto impl(urng_t && urange)
    {
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The views::all adaptor can only be passed viewable_ranges, i.e. Views or &-to-non-View.");

        // views are always passed as-is
        if constexpr (std::ranges::view<remove_cvref_t<urng_t>>)
        {
            return std::views::all(std::forward<urng_t>(urange));
        }
        // string const &
        else if constexpr (is_type_specialisation_of_v<remove_cvref_t<urng_t>, std::basic_string> &&
                           std::is_const_v<std::remove_reference_t<urng_t>>)
        {
            return std::basic_string_view{std::ranges::data(urange), std::ranges::size(urange)};
        }
        // contiguous
        else if constexpr (forwarding_range<urng_t> &&
                           std::ranges::contiguous_range<urng_t> &&
                           std::ranges::sized_range<urng_t>)
        {
            return std::span{std::ranges::data(urange), std::ranges::size(urange)};
        }
        // random_access
        else if constexpr (forwarding_range<urng_t> &&
                           std::ranges::random_access_range<urng_t> &&
                           std::ranges::sized_range<urng_t>)
        {
            return std::ranges::subrange<std::ranges::iterator_t<urng_t>, std::ranges::iterator_t<urng_t>>
            {
                std::ranges::begin(urange),
                std::ranges::begin(urange) + std::ranges::size(urange),
                std::ranges::size(urange)
            };
        }
        // pass to std::views::all (will return ref-view)
        else
        {
            return std::views::all(std::forward<urng_t>(urange));
        }
    }
};

} // namespace seqan3::detail

// ============================================================================
//  views::all (adaptor instance definition)
// ============================================================================

namespace seqan3::views
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that behaves like std::views::all, but type erases contiguous ranges.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             Either `std::basic_string_view{urange}`, or `std::span{urange}`, or `std::views::all(urange)`.
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/view_all.hpp}
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)  | `rrng_t` (returned range type)  |
 * |----------------------------------|:---------------------------------:|:-------------------------------:|
 * | std::ranges::input_range         | *required*                        | *preserved*                     |
 * | std::ranges::forward_range       |                                   | *preserved*                     |
 * | std::ranges::bidirectional_range |                                   | *preserved*                     |
 * | std::ranges::random_access_range |                                   | *preserved*                     |
 * | std::ranges::contiguous_range    |                                   | *preserved*                     |
 * |                                  |                                   |                                  |
 * | std::ranges::viewable_range      | *required*                        | *guaranteed*                    |
 * | std::ranges::view                |                                   | *guaranteed*                    |
 * | std::ranges::sized_range         |                                   | *preserved*                     |
 * | std::ranges::common_range        |                                   | *preserved*                     |
 * | std::ranges::output_range        |                                   | *preserved*                     |
 * | seqan3::const_iterable_range     |                                   | *preserved*                     |
 * |                                  |                                   |                                  |
 * | std::ranges::range_reference_t   |                                   | seqan3::reference_t<urng_t>     |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Return type
 *
 * | `urng_t` (underlying range type)                                                       | `rrng_t` (returned range type)  |
 * |:--------------------------------------------------------------------------------------:|:-------------------------------:|
 * | `std::basic_string const &` *or* `std::basic_string_view`                              | `std::basic_string_view`        |
 * | `seqan3::forwarding_range && std::ranges::sized_range && std::ranges::contiguous_range`   | `std::span`                     |
 * | `seqan3::forwarding_range && std::ranges::sized_range && std::ranges::random_access_range` | `std::ranges::subrange`         |
 * | *else*                                                                                 | *implementation defined type*   |
 *
 * This adaptor is different from std::views::take in that it performs type erasure for some underlying ranges.
 * It returns exactly the type specified above.
 *
 * ### Example
 *
* \include test/snippet/range/views/view_all.cpp
 *
 * \hideinitializer
 */
inline constexpr auto all = detail::all_fn{};

//!\}

} // namespace seqan3::views

namespace seqan3
{
//!\brief Deduces the return value of seqan3::views::all.
template <typename t>
using all_view = decltype(views::all(std::declval<t>()));
}
