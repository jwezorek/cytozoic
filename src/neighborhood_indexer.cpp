#include "neighborhood_indexer.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <map>
#include <vector>
#include <ranges>

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
        )
    };

}

cz::sum_of_states_indexer::sum_of_states_indexer(std::size_t max_neighbors)
    : max_neighbors_(max_neighbors)
{}

std::size_t cz::sum_of_states_indexer::column_index( const std::vector<int8_t>& neighbor_states,
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
        neighbor_states.end());

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

std::size_t cz::max_state_indexer::num_columns( std::size_t num_states) const {
    if (num_states == 0) {
        throw std::runtime_error(
            "max_state_indexer requires num_states > 0."
        );
    }

    return num_states;
}

std::string cz::max_state_indexer::name() const {
    return "max state";
}

std::vector<std::string> cz::named_indexers() {
    return g_indexers | rv::transform(
        [](const auto& indexer) {
            return indexer->name();
        }
    ) | r::to<std::vector>();
}

cz::neighborhood_indexer cz::indexer_from_name(const std::string& str) {
    return{};
}
