#pragma once

#include "cytozoic.hpp"
#include <string>
#include <optional>

/*------------------------------------------------------------------------------------------------*/

namespace cz {

    bool save_ruleset_to_file(
        const std::string& file_path,
        const cyto_params& params
    );

    std::optional<cyto_params> load_ruleset_from_file(const std::string& file_path);

}