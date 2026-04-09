#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "third-party/proxy/proxy.h"

/*------------------------------------------------------------------------------------------------*/

namespace cz {

    PRO_DEF_MEM_DISPATCH(mem_column_index, column_index);
    PRO_DEF_MEM_DISPATCH(mem_num_columns, num_columns);

    struct neighborhood_indexer_facade : pro::facade_builder
        ::add_convention<
            mem_column_index,
            std::size_t(
                const std::vector<int8_t>& neighbor_states,
                std::size_t num_states
            ) const
        >
        ::add_convention<
            mem_num_columns,
            std::size_t(std::size_t num_states) const
        >
        ::support_copy<pro::constraint_level::nontrivial>
        ::build {};

    using neighborhood_indexer = pro::proxy<neighborhood_indexer_facade>;

    struct sum_of_states_indexer {
        explicit sum_of_states_indexer(std::size_t max_neighbors);

        std::size_t column_index(
            const std::vector<int8_t>& neighbor_states,
            std::size_t num_states) const;

        std::size_t num_columns(std::size_t num_states) const;

    private:
        std::size_t max_neighbors_;
    };

    struct max_state_indexer {
        std::size_t column_index(
            const std::vector<int8_t>& neighbor_states,
            std::size_t num_states) const;

        std::size_t num_columns(std::size_t num_states) const;
    };

} // namespace cz