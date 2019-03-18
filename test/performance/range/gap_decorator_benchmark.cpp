// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cstdlib>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <utility>

#include <benchmark/benchmark.h>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/decorator/all.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{
/* Mocked gap decorator using a modifyable container of the union type of gap and
 * alphabet. This class serves for comparison with decorators that take a reference
 * to the underlying sequence and do not modify it.
 */
template<std::ranges::ViewableRange inner_type>
requires std::ranges::RandomAccessRange<inner_type> && std::ranges::SizedRange<inner_type> &&
         (std::is_const_v<std::remove_reference_t<inner_type>> || std::ranges::View<inner_type>)
class gapped_sequence
{
private:
    using inner_sequence_type = typename std::vector<gapped<value_type_t<inner_type>>>;
    inner_sequence_type gapseq;

    class gapped_sequence_iterator
    {
    private:
        typename std::add_pointer_t<gapped_sequence const> host{nullptr};
        void jump(typename gapped_sequence::size_type const new_pos)
        {
            assert(new_pos <= host->size());
            pos = new_pos;
        }
    public:
        typename gapped_sequence::size_type pos{0u};

        using difference_type = typename gapped_sequence::difference_type;
        using value_type = typename gapped_sequence::value_type;
        using reference = typename gapped_sequence::const_reference;
        using const_reference = reference;
        using pointer = value_type *;
        using iterator_category = std::bidirectional_iterator_tag;

        constexpr gapped_sequence_iterator() = default;
        constexpr gapped_sequence_iterator(gapped_sequence_iterator const &) = default;
        constexpr gapped_sequence_iterator & operator=(gapped_sequence_iterator const &) = default;
        constexpr gapped_sequence_iterator (gapped_sequence_iterator &&) = default;
        constexpr gapped_sequence_iterator & operator=(gapped_sequence_iterator &&) = default;
        ~gapped_sequence_iterator() = default;

        explicit constexpr gapped_sequence_iterator(gapped_sequence const & host_) : host(&host_) {}

        constexpr gapped_sequence_iterator(gapped_sequence const & host_,
            typename gapped_sequence::size_type const pos_) : host(&host_)
        {
            jump(pos_);
        }

        constexpr gapped_sequence_iterator & operator++() noexcept
        {
            ++pos;
            return *this;
        }

        constexpr gapped_sequence_iterator & operator--() noexcept
        {
            --pos;
            return *this;
        }

        constexpr gapped_sequence_iterator operator++(int) noexcept
        {
            gapped_sequence_iterator cpy{*this};
            ++(*this);
            return cpy;
        }

        constexpr gapped_sequence_iterator operator--(int) noexcept
        {
            gapped_sequence_iterator cpy{*this};
            --(*this);
            return cpy;
        }

        constexpr reference operator*() const noexcept
        {
            return host->at(pos);
        }

        constexpr pointer operator->() const noexcept
        {
            return &(host->at(pos));
        }

        constexpr friend bool operator==(gapped_sequence_iterator const & lhs,
                                         gapped_sequence_iterator const & rhs)
        {
            return lhs.pos == rhs.pos;
        }

        constexpr friend bool operator!=(gapped_sequence_iterator const & lhs,
                                         gapped_sequence_iterator const & rhs)
        {
            return lhs.pos != rhs.pos;
        }

        constexpr friend bool operator<(gapped_sequence_iterator const & lhs,
                                        gapped_sequence_iterator const & rhs)
        {
            return lhs.pos < rhs.pos;
        }

        constexpr friend bool operator>(gapped_sequence_iterator const & lhs,
                                        gapped_sequence_iterator const & rhs)
        {
            return lhs.pos > rhs.pos;
        }

        constexpr friend bool operator<=(gapped_sequence_iterator const & lhs,
                                         gapped_sequence_iterator const & rhs)
        {
            return lhs.pos <= rhs.pos;
        }

        constexpr friend bool operator>=(gapped_sequence_iterator const & lhs,
                                         gapped_sequence_iterator const & rhs)
        {
            return lhs.pos >= rhs.pos;
        }
    };
public:
    using value_type = gapped<value_type_t<inner_type>>;
    using reference = value_type;
    using const_reference = reference;
    using size_type = size_type_t<inner_type>;
    using difference_type = difference_type_t<inner_type>;
    using iterator = gapped_sequence_iterator;
    using const_iterator = iterator;
    using ungapped_range_type = inner_type;

    constexpr gapped_sequence() = default;
    constexpr gapped_sequence(gapped_sequence const &) = default;
    constexpr gapped_sequence & operator=(gapped_sequence const &) = default;
    constexpr gapped_sequence(gapped_sequence && rhs) = default;
    constexpr gapped_sequence & operator=(gapped_sequence && rhs) = default;
    ~gapped_sequence() = default;
    // upon construction by reference create local sequence
    constexpr gapped_sequence(inner_type const & range)
    {
        for (auto elem : range)
            gapseq.push_back(gapped<value_type_t<inner_type>>{elem});
    }

    size_type size() const noexcept
    {
        return gapseq.size();
    }

    iterator insert_gap(iterator const it, size_type const count = 1)
    {
        std::vector<gap> gaps(count, gap{});
        auto gapseq_it = gapseq.begin() + it.pos;
        gapseq.insert(gapseq_it, gaps.begin(), gaps.end());
        return it;
    }

    iterator erase_gap(iterator const it)
    {
        if ((*it) != gap{})
            throw gap_erase_failure{"The range to be erased does not correspond to a consecutive gap."};
        auto end_it = std::next(it);
        return erase_gap(it, end_it);
    }

    iterator erase_gap(iterator const first, iterator const last)
    {
        for (auto it = first; it < last; ++it)
        {
            if (*it != gap{})
                throw gap_erase_failure{"There is no gap to erase in range."};
        }
        typename inner_sequence_type::iterator it1 = gapseq.begin() + first.pos;
        typename inner_sequence_type::iterator it2 = gapseq.begin() + last.pos;
        gapseq.erase(it1, it2);
        return first;
    }

    iterator begin() const noexcept
    {
        return iterator{*this};
    }

    const_iterator cbegin() const noexcept
    {
        return const_iterator{*this};
    }

    reference at(size_type const i)
    {
        if (i >= size())
            throw std::out_of_range{"Trying to access element behind the last in gap_decorator."};
        return (*this)[i];
    }

    const_reference at(size_type const i) const
    {
        if (i >= size())
            throw std::out_of_range{"Trying to access element behind the last in gap_decorator."};
        return (*this)[i];
    }

    constexpr reference operator[](size_type const i) const noexcept
    {
        return gapseq[i];
    }
};

} // namespace seqan3

using namespace seqan3;

typedef long double time_type;
/* Helper function to sample gap length for each ungapped sequence position.
 *
 * Parameter:
 * size_type    size type of underlying sequence range
 *
 * Arguments:
 * gap_vector   reference to empty vector storing the sampled gap lengths
 * size         final (aligned) sequence size
 */
template<typename size_type>
void sample(std::vector<size_type> * gap_vector, size_type size)
{
    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::array<double,10> cumsum = {0.6395, 0.8263, 0.8871, 0.9257, 0.9544, 0.9709, 0.9813, 0.9890, 0.9955, 1.0000};
    for (size_type i = 0; i < size; ++i){
        double y = uni(generator);
        gap_vector->at(i) = y;
        auto it = std::find_if(cumsum.begin(), cumsum.end(), [y](double i){return y <= i;});
        gap_vector->at(i) = it - cumsum.begin();
    }
}

/* Helper function to adjust the ungapped sequence length w.r.t. sampled gaps s.t.
 * the gapped sequence length does not exceed the targeted length.
 *
 * Parameters:
 * size_type        size type of the sequence
 * sequence_type    ungapped sequence type
 *
 * Arguments:
 * gaps         reference to the gap vector
 * seq          reference to ungapped sequence
 * seq_len      final sequence length
 */
template<typename size_type, typename sequence_type>
void resize(std::vector<size_type> & gaps, sequence_type & seq, unsigned int seq_len)
{
    size_type letter_acc = 0;
    size_type gap_pos = 0;
    size_type gap_acc = 0;

    while (gap_pos < gaps.size() && gap_acc + letter_acc < seq_len)
    {
        if (!gaps[gap_pos])
            ++letter_acc;
        else
        {
            if (letter_acc + gap_acc + gaps[gap_pos] > seq_len)
            {
                gaps[gap_pos] = seq_len - gap_acc - letter_acc;
                gap_acc += gaps[gap_pos];
                ++gap_pos;
                break;
            }
            else
                gap_acc += gaps[gap_pos];
        }
        ++gap_pos;
    }
    seq.resize(std::max<size_type>(1, letter_acc));  // resize ungapped sequence
    gaps.resize(gap_pos);      // trim sampled gap vector
}

/* Helper function to prepare a gapped sequence for the benchmark (case gap_flag=true)
 *
 * Parameters:
 * size_type        size type of the sequence
 * iterator_type    iterator type of decorated sequence iterator
 * gap_decorator_t  gap decorator type, e.g. gap_decorator_anchor_set
 *
 * Arguments:
 * gaps             reference to gap vector
 * gap_decorator    reference to gap decorator
 */
template<typename size_type, typename iterator_type, typename gap_decorator_t>
void insert_gaps(std::vector<size_type> & gaps, gap_decorator_t & gap_decorator)
{
    [[maybe_unused]] size_type gap_acc = 0;
    for (size_type i = 0; i < gaps.size(); ++i)
    {
        if (gaps[i])
        {
            iterator_type it(gap_decorator, std::min(i + gap_acc, gap_decorator.size()));
            gap_decorator.insert_gap(it, gaps[i]);
        }
        gap_acc += gaps[i];
    }
}

// ============================================================================
//  read left to right (looped in case #ops exceeds sequence length)
// ============================================================================
/* Parameters:
 * gap_decorator_t      gap decorator class, e.g. gap_decorator_anchor_set
 * ungapped_range_t     unaligned sequence to which a reference is stored
 * alphabet_t           alphabet type of the ungapped sequence
 * seq_len              final sequence length including gap symbols
 * gapped_flag          operate on already gapped (true) or ungapped sequence (false)
 */
template <template <typename> typename gap_decorator_t,
       template<typename> typename ungapped_range_t, typename alphabet_t, bool gapped_flag>
static void read_left2right(benchmark::State& state)
{
    using range_reference_t = const ungapped_range_t<alphabet_t> &;
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t<range_reference_t>::size_type;
    using iterator_type = typename gap_decorator_t<range_reference_t>::iterator;
    using sequence_type = ungapped_range_t<alphabet_t>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(&gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t<range_reference_t> gap_decorator(seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, iterator_type, gap_decorator_t<range_reference_t>>(gaps, gap_decorator);

    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pos = op_ctr % seq_len;
        state.ResumeTiming();
        benchmark::DoNotOptimize(gap_decorator[pos]);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["read_op"] = op_ctr;
}

// 1 a) Read from left to right in ungapped sequence
BENCHMARK_TEMPLATE(read_left2right, gap_decorator_anchor_set, std::vector, dna4, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(read_left2right, gapped_sequence, std::vector, dna4, false)->Range(1<<2, 1<<15);
// 1 b) Read from left to right in gapped sequence
BENCHMARK_TEMPLATE(read_left2right, gap_decorator_anchor_set, std::vector, dna4, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(read_left2right, gapped_sequence, std::vector, dna4, true)->Range(1<<2, 1<<15);


// ============================================================================
//  read at random position
// ============================================================================
/* Parameters:
 * gap_decorator_t      gap decorator class, e.g. gap_decorator_anchor_set
 * ungapped_range_t     unaligned sequence to which a reference is stored
 * alphabet_t           alphabet type of the ungapped sequence
 * seq_len              final sequence length including gap symbols
 * gapped_flag          operate on already gapped (true) or ungapped sequence (false)
 */

template <template <typename> typename gap_decorator_t,
          template<typename> typename ungapped_range_t, typename alphabet_t, bool gapped_flag>
static void read_random(benchmark::State& state)
{
    using range_reference_t = const ungapped_range_t<alphabet_t> &;
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t<range_reference_t>::size_type;
    using iterator_type = typename gap_decorator_t<range_reference_t>::iterator;
    using sequence_type = ungapped_range_t<alphabet_t>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(&gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t<range_reference_t> gap_decorator(seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, iterator_type, gap_decorator_t<range_reference_t>>(gaps, gap_decorator);

    std::mt19937 generator(time(0)); //Standard mersenne_twister_engine seeded with current time
    std::uniform_real_distribution<> uni_dis(0.0, static_cast<double>(seq_len));

    size_t pos = 0;
    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        pos = uni_dis(generator);
        state.ResumeTiming();
        benchmark::DoNotOptimize(gap_decorator[pos]);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["read_op"] = op_ctr;
}


// 2 a) Read at random position in ungapped sequence
BENCHMARK_TEMPLATE(read_random, gap_decorator_anchor_set, std::vector, dna4, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(read_random, gapped_sequence, std::vector, dna4, false)->Range(1<<2, 1<<15);
// 2 b) Read at random position in gapped sequence
BENCHMARK_TEMPLATE(read_random, gap_decorator_anchor_set, std::vector, dna4, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(read_random, gapped_sequence, std::vector, dna4, true)->Range(1<<2, 1<<15);

// ============================================================================
//  insert left to right
// ============================================================================
template <template <typename> typename gap_decorator_t,
       template<typename> typename ungapped_range_t, typename alphabet_t, bool gapped_flag>
static void insert_left2right(benchmark::State& state)
{
    using range_reference_t = const ungapped_range_t<alphabet_t> &;
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t<range_reference_t>::size_type;
    using iterator_type = typename gap_decorator_t<range_reference_t>::iterator;
    using sequence_type = ungapped_range_t<alphabet_t>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(&gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t<range_reference_t> gap_decorator(seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, iterator_type, gap_decorator_t<range_reference_t>>(gaps, gap_decorator);

    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pos = op_ctr % seq_len;
        iterator_type it(gap_decorator, pos);
        state.ResumeTiming();
        gap_decorator.insert_gap(it, 1);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["insert_op"] = op_ctr;
}

// 3 a) Insert gaps of length 1 from left to right into ungapped sequence
BENCHMARK_TEMPLATE(insert_left2right, gap_decorator_anchor_set, std::vector, dna4, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_left2right, gapped_sequence, std::vector, dna4, false)->Range(1<<2, 1<<15);
// 3 b) Insert gaps of length 1 from left to right into gapped sequence
BENCHMARK_TEMPLATE(insert_left2right, gap_decorator_anchor_set, std::vector, dna4, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_left2right, gapped_sequence, std::vector, dna4, true)->Range(1<<2, 1<<15);

// ============================================================================
//  insert right to left
// ============================================================================
template <template <typename> typename gap_decorator_t,
       template<typename> typename ungapped_range_t, typename alphabet_t, bool gapped_flag>
static void insert_right2left(benchmark::State& state)
{
    using range_reference_t = const ungapped_range_t<alphabet_t> &;
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t<range_reference_t>::size_type;
    using iterator_type = typename gap_decorator_t<range_reference_t>::iterator;
    using sequence_type = ungapped_range_t<alphabet_t>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(&gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t<range_reference_t> gap_decorator(seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, iterator_type, gap_decorator_t<range_reference_t>>(gaps, gap_decorator);

    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pos = seq_len - (op_ctr % seq_len) - 1;
        iterator_type it(gap_decorator, pos);
        state.ResumeTiming();
        gap_decorator.insert_gap(it, 1);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["insert_op"] = op_ctr;
}

// 4 a) Insert gaps of length 1 from left to right into ungapped sequence
BENCHMARK_TEMPLATE(insert_right2left, gap_decorator_anchor_set, std::vector, dna4, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_right2left, gapped_sequence, std::vector, dna4, false)->Range(1<<2, 1<<15);
// 4 b) Insert gaps of length 1 from left to right into gapped sequence
BENCHMARK_TEMPLATE(insert_right2left, gap_decorator_anchor_set, std::vector, dna4, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_right2left, gapped_sequence, std::vector, dna4, true)->Range(1<<2, 1<<15);

// ============================================================================
//  insert at random position
// ============================================================================
template <template <typename> typename gap_decorator_t,
       template<typename> typename ungapped_range_t, typename alphabet_t, bool gapped_flag>
static void insert_random(benchmark::State& state)
{
    using range_reference_t = const ungapped_range_t<alphabet_t> &;
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t<range_reference_t>::size_type;
    using iterator_type = typename gap_decorator_t<range_reference_t>::iterator;
    using sequence_type = ungapped_range_t<alphabet_t>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(&gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t<range_reference_t> gap_decorator(seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, iterator_type, gap_decorator_t<range_reference_t>>(gaps, gap_decorator);

    std::mt19937 generator(time(0)); //Standard mersenne_twister_engine seeded with current time
    std::uniform_real_distribution<> uni_dis(0.0, static_cast<double>(seq_len));

    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pos = uni_dis(generator);
        iterator_type it(gap_decorator, pos);
        state.ResumeTiming();
        gap_decorator.insert_gap(it, 1);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["insert_op"] = op_ctr;
}

// 5 a) Insert gaps of length 1 at random position into ungapped sequence
BENCHMARK_TEMPLATE(insert_random, gap_decorator_anchor_set, std::vector, dna4, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_random, gapped_sequence, std::vector, dna4, false)->Range(1<<2, 1<<15);
// 5 b) Insert gaps of length 1 at random position into gapped sequence
BENCHMARK_TEMPLATE(insert_random, gap_decorator_anchor_set, std::vector, dna4, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_random, gapped_sequence, std::vector, dna4, true)->Range(1<<2, 1<<15);

// ============================================================================
//  delete at random position
// ============================================================================
template <template <typename> typename gap_decorator_t,
       template<typename> typename ungapped_range_t, typename alphabet_t, bool gapped_flag>
static void delete_random(benchmark::State& state)
{
    using range_reference_t = const ungapped_range_t<alphabet_t> &;
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t<range_reference_t>::size_type;
    using iterator_type = typename gap_decorator_t<range_reference_t>::iterator;
    using sequence_type = ungapped_range_t<alphabet_t>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(&gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t<range_reference_t> gap_decorator(seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, iterator_type, gap_decorator_t<range_reference_t>>(gaps, gap_decorator);

    std::mt19937 generator(time(0)); //Standard mersenne_twister_engine seeded with current time
    std::uniform_real_distribution<> uni_dis(0.0, static_cast<double>(seq_len));

    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pos = uni_dis(generator);
        iterator_type first(gap_decorator, pos);
        gap_decorator.insert_gap(first, 2);
        iterator_type last(gap_decorator, pos+2);
        state.ResumeTiming();
        gap_decorator.erase_gap(first, last);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["delete_op"] = op_ctr;
}

// 6 a) Erase gaps at random position from initially ungapped sequence
BENCHMARK_TEMPLATE(delete_random, gap_decorator_anchor_set, std::vector, dna4, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(delete_random, gapped_sequence, std::vector, dna4, false)->Range(1<<2, 1<<15);
// 6 b) Erase gaps at random position from initially gapped sequence
BENCHMARK_TEMPLATE(delete_random, gap_decorator_anchor_set, std::vector, dna4, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(delete_random, gapped_sequence, std::vector, dna4, true)->Range(1<<2, 1<<15);

BENCHMARK_MAIN();