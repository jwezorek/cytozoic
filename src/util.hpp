#pragma once

#include "types.hpp"
#include <vector>
#include <optional>

namespace cz {

    std::vector<point> random_points(size_t n,
        double bounds_wd = 1.0, double bounds_hgt = 1.0,
        std::optional<uint64_t> seed = {}
    );

}