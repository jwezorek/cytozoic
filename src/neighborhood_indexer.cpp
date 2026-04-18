#include "neighborhood_indexer.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {

    std::vector<cz::neighborhood_indexer> g_indexers = {

        pro::make_proxy<cz::neighborhood_indexer_facade>(
            cz::sum_of_states_indexer{10}
        ),

        pro::make_proxy<cz::neighborhood_indexer_facade>(
            cz::max_state_indexer{}
        ),

        pro::make_proxy<cz::neighborhood_indexer_facade>(
            cz::trinary_histogram_indexer{}
        ),

        pro::make_proxy<cz::neighborhood_indexer_facade>(
            cz::quarternary_histogram_indexer{}
        ),

        pro::make_proxy<cz::neighborhood_indexer_facade>(
            cz::min_max_state_indexer{}
        )
    };

    int to_trinary_digit(int count) {
        if (count <= 1) {
            return 0;
        }
        else if (count <= 3) {
            return 1;
        }

        return 2;
    }

    int to_quarternary_digit(int count) {
        if (count <= 0) {
            return 0;
        }
        else if (count == 1) {
            return 1;
        }
        else if (count <= 3) {
            return 2;
        }

        return 3;
    }

    std::size_t a_to_the_nth(int a, std::size_t n) {
        std::size_t value = 1;

        for (std::size_t i = 0; i < n; ++i) {
            value *= a;
        }

        return value;
    }

    std::size_t trinary_to_integer(const std::vector<int>& digits) {
        std::size_t value = 0;
        std::size_t place = 1;

        for (int digit : digits) {
            value += static_cast<std::size_t>(digit) * place;
            place *= 3;
        }

        return value;
    }

    std::size_t quarternary_to_integer(const std::vector<int>& digits) {
        std::size_t value = 0;
        std::size_t place = 1;

        for (int digit : digits) {
            value += static_cast<std::size_t>(digit) * place;
            place *= 4;
        }

        return value;
    }

} // namespace

/*------------------------------------------------------------------------------------------------*/

cz::sum_of_states_indexer::sum_of_states_indexer(std::size_t max_neighbors)
    : max_neighbors_(max_neighbors)
{}

std::size_t cz::sum_of_states_indexer::column_index(
    const std::vector<int8_t>& neighbor_states,
    std::size_t num_states) const
{
    if (num_states == 0) {
        throw std::runtime_error(
            "sum_of_states_indexer requires num_states > 0."
        );
    }

    const std::size_t max_sum = max_neighbors_ * (num_states - 1);

    const std::size_t sum = std::transform_reduce(
        neighbor_states.begin(),
        neighbor_states.end(),
        std::size_t{ 0 },
        std::plus<>{},
        [num_states](int8_t state) -> std::size_t {
            if (state < 0) {
                throw std::runtime_error(
                    "sum_of_states_indexer requires non-negative states."
                );
            }

            const std::size_t value = static_cast<std::size_t>(state);
            if (value >= num_states) {
                throw std::runtime_error(
                    "sum_of_states_indexer received a state outside "
                    "the valid range."
                );
            }

            return value;
        }
    );

    return std::min(sum, max_sum);
}

std::size_t cz::sum_of_states_indexer::num_columns(
    std::size_t num_states) const
{
    if (num_states == 0) {
        throw std::runtime_error(
            "sum_of_states_indexer requires num_states > 0."
        );
    }

    return max_neighbors_ * (num_states - 1) + 1;
}

std::string cz::sum_of_states_indexer::name() const
{
    return "sum of states";
}

/*------------------------------------------------------------------------------------------------*/

std::size_t cz::max_state_indexer::column_index(
    const std::vector<int8_t>& neighbor_states,
    std::size_t num_states) const
{
    if (neighbor_states.empty()) {
        throw std::runtime_error(
            "max_state_indexer requires a non-empty state vector."
        );
    }

    if (num_states == 0) {
        throw std::runtime_error(
            "max_state_indexer requires num_states > 0."
        );
    }

    const auto it = std::max_element(
        neighbor_states.begin(),
        neighbor_states.end()
    );

    if (*it < 0) {
        throw std::runtime_error(
            "max_state_indexer requires non-negative states."
        );
    }

    const std::size_t value = static_cast<std::size_t>(*it);
    if (value >= num_states) {
        throw std::runtime_error(
            "max_state_indexer received a state outside the valid range."
        );
    }

    return value;
}

std::size_t cz::max_state_indexer::num_columns(std::size_t num_states) const
{
    if (num_states == 0) {
        throw std::runtime_error(
            "max_state_indexer requires num_states > 0."
        );
    }

    return num_states;
}

std::string cz::max_state_indexer::name() const
{
    return "max state";
}

/*------------------------------------------------------------------------------------------------*/

std::size_t cz::trinary_histogram_indexer::column_index(
    const std::vector<int8_t>& neighbor_states,
    std::size_t num_states) const
{
    if (num_states == 0) {
        throw std::runtime_error(
            "trinary_histogram_indexer requires num_states > 0."
        );
    }

    std::vector<int> histogram(num_states, 0);

    for (int8_t s : neighbor_states) {
        if (s < 0) {
            throw std::runtime_error(
                "trinary_histogram_indexer requires non-negative states."
            );
        }

        const auto state = static_cast<std::size_t>(s);

        if (state >= num_states) {
            throw std::runtime_error(
                "trinary_histogram_indexer received a state outside the valid range."
            );
        }

        ++histogram[state];
    }

    const auto digits =
        histogram
        | rv::transform(to_trinary_digit)
        | r::to<std::vector>();

    return trinary_to_integer(digits);
}

std::size_t cz::trinary_histogram_indexer::num_columns(
    std::size_t num_states) const
{
    if (num_states == 0) {
        throw std::runtime_error(
            "trinary_histogram_indexer requires num_states > 0."
        );
    }

    return a_to_the_nth(3, num_states);
}

std::string cz::trinary_histogram_indexer::name() const
{
    return "trinary histogram";
}

/*------------------------------------------------------------------------------------------------*/

std::size_t cz::quarternary_histogram_indexer::column_index(
    const std::vector<int8_t>& neighbor_states,
    std::size_t num_states) const
{
    if (num_states == 0) {
        throw std::runtime_error(
            "quarternary_histogram_indexer requires num_states > 0."
        );
    }

    std::vector<int> histogram(num_states, 0);

    for (int8_t s : neighbor_states) {
        if (s < 0) {
            throw std::runtime_error(
                "quarternary_histogram_indexer requires non-negative states."
            );
        }

        const auto state = static_cast<std::size_t>(s);

        if (state >= num_states) {
            throw std::runtime_error(
                "quarternary_histogram_indexer received a state outside the valid range."
            );
        }

        ++histogram[state];
    }

    const auto digits =
        histogram | rv::transform(to_quarternary_digit) | r::to<std::vector>();

    return quarternary_to_integer(digits);
}

std::size_t cz::quarternary_histogram_indexer::num_columns(
    std::size_t num_states) const
{
    if (num_states == 0) {
        throw std::runtime_error(
            "quarternary_histogram_indexer requires num_states > 0."
        );
    }

    return a_to_the_nth(4, num_states);
}

std::string cz::quarternary_histogram_indexer::name() const
{
    return "quarternary histogram";
}

/*------------------------------------------------------------------------------------------------*/

std::vector<std::string> cz::named_indexers()
{
    return g_indexers | rv::transform(
        [](const auto& indexer) {
            return indexer->name();
        }
    ) | r::to<std::vector>();
}

cz::neighborhood_indexer cz::indexer_from_name(const std::string& str)
{
    static std::map<std::string, cz::neighborhood_indexer> name_to_indexer;

    if (name_to_indexer.empty()) {
        for (const auto& indexer : g_indexers) {
            name_to_indexer[indexer->name()] = indexer;
        }
    }

    return name_to_indexer.at(str);
}

std::size_t cz::min_max_state_indexer::column_index(
    const std::vector<int8_t>& neighbor_states,
    std::size_t num_states) const {
    if (neighbor_states.empty()) {
        throw std::runtime_error(
            "min_max_state_indexer requires a non-empty state vector."
        );
    }

    if (num_states == 0) {
        throw std::runtime_error(
            "min_max_state_indexer requires num_states > 0."
        );
    }

    const auto [min_it, max_it] = std::minmax_element(
        neighbor_states.begin(),
        neighbor_states.end()
    );

    if (*min_it < 0 || *max_it < 0) {
        throw std::runtime_error(
            "min_max_state_indexer requires non-negative states."
        );
    }

    const std::size_t min_state = static_cast<std::size_t>(*min_it);
    const std::size_t max_state = static_cast<std::size_t>(*max_it);

    if (max_state >= num_states) {
        throw std::runtime_error(
            "min_max_state_indexer received a state outside the valid range."
        );
    }

    if (min_state > max_state) {
        throw std::runtime_error(
            "min_max_state_indexer computed invalid min/max ordering."
        );
    }

    const std::size_t prefix =
        min_state * num_states - (min_state * (min_state - 1)) / 2;

    return prefix + (max_state - min_state);
}

std::size_t cz::min_max_state_indexer::num_columns(
    std::size_t num_states) const {
    if (num_states == 0) {
        throw std::runtime_error(
            "min_max_state_indexer requires num_states > 0."
        );
    }

    return num_states * (num_states + 1) / 2;
}

std::string cz::min_max_state_indexer::name() const
{
    return "min/max state";
}
