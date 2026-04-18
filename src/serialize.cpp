#include "serialize.hpp"
#include "third-party/json.hpp"
#include <fstream>

/*------------------------------------------------------------------------------------------------*/

namespace {

    using json = nlohmann::json;

    json color_to_json(const cz::color& color)
    {
        return json{
            { "r", color.r },
            { "g", color.g },
            { "b", color.b }
        };
    }

    cz::color color_from_json(const json& j)
    {
        return cz::color{
            static_cast<uint8_t>(j.at("r").get<int>()),
            static_cast<uint8_t>(j.at("g").get<int>()),
            static_cast<uint8_t>(j.at("b").get<int>())
        };
    }

    json state_table_row_to_json(const cz::state_table_row& row)
    {
        json j = json::array();

        for (int8_t value : row) {
            j.push_back(static_cast<int>(value));
        }

        return j;
    }

    cz::state_table_row state_table_row_from_json(const json& j)
    {
        if (!j.is_array()) {
            throw std::runtime_error("state table row must be a JSON array.");
        }

        cz::state_table_row row;
        row.reserve(j.size());

        for (const auto& value : j) {
            row.push_back(static_cast<int8_t>(value.get<int>()));
        }

        return row;
    }

    json state_table_to_json(const cz::state_table& table)
    {
        json j = json::array();

        for (const auto& row : table) {
            j.push_back(state_table_row_to_json(row));
        }

        return j;
    }

    cz::state_table state_table_from_json(const json& j)
    {
        if (!j.is_array()) {
            throw std::runtime_error("state table must be a JSON array.");
        }

        cz::state_table table;
        table.reserve(j.size());

        for (const auto& row : j) {
            table.push_back(state_table_row_from_json(row));
        }

        return table;
    }

    std::vector<double> double_vector_from_json(const json& j)
    {
        if (!j.is_array()) {
            throw std::runtime_error("expected a JSON array of numbers.");
        }

        std::vector<double> values;
        values.reserve(j.size());

        for (const auto& value : j) {
            values.push_back(value.get<double>());
        }

        return values;
    }

    json double_vector_to_json(const std::vector<double>& values)
    {
        json j = json::array();

        for (double value : values) {
            j.push_back(value);
        }

        return j;
    }

    cz::color_table color_table_from_json(const json& j)
    {
        if (!j.is_array()) {
            throw std::runtime_error("palette must be a JSON array.");
        }

        cz::color_table palette;
        palette.reserve(j.size());

        for (const auto& value : j) {
            palette.push_back(color_from_json(value));
        }

        return palette;
    }

    json color_table_to_json(const cz::color_table& palette)
    {
        json j = json::array();

        for (const auto& color : palette) {
            j.push_back(color_to_json(color));
        }

        return j;
    }

}

bool cz::save_ruleset_to_file( const std::string& file_path, const cyto_params& params ) {

    try {

        json j = {
            { "format", "cytozoic_ruleset" },
            { "version", 1 },
            { "num_states", params.num_states },
            { "num_initial_cells", params.num_initial_cells },
            { "initial_state_density", double_vector_to_json(params.initial_state_density) },
            { "palette", color_table_to_json(params.palette) },
            { "cell_indexer", params.cell_indexer->name() },
            { "cell_state_table", state_table_to_json(params.cell_state_table) },
            { "vertex_indexer", params.vertex_indexer->name() },
            { "vertex_table", state_table_row_to_json(params.vertex_table) }
        };

        std::ofstream file(file_path, std::ios::binary | std::ios::trunc);
        if (!file) {
            return false;
        }

        file << j.dump(4);
        file.close();

        return true;

    } catch (const std::exception& ex) {

        return false;
    }

}

std::optional<cz::cyto_params> cz::load_ruleset_from_file(const std::string& file_path) {

    try {
        std::ifstream file(file_path, std::ios::binary);

        if (!file) {
            return {};
        }

        const json j = json::parse(file);

        if (!j.is_object()) {
            throw std::runtime_error(
                "ruleset file must contain a top-level JSON object."
            );
        }

        const std::string format = j.value("format", std::string{});
        const int version = j.value("version", 0);

        if (format != "cytozoic_ruleset") {
            throw std::runtime_error("unrecognized ruleset format.");
        }

        if (version != 1) {
            throw std::runtime_error("unsupported ruleset version.");
        }

        cz::cyto_params loaded;

        loaded.num_states = j.at("num_states").get<int>();
        loaded.num_initial_cells = j.at("num_initial_cells").get<int>();
        loaded.initial_state_density =
            double_vector_from_json(j.at("initial_state_density"));
        loaded.palette =
            color_table_from_json(j.at("palette"));

        loaded.cell_indexer = cz::indexer_from_name(
            j.at("cell_indexer").get<std::string>()
        );

        loaded.vertex_indexer = cz::indexer_from_name(
            j.at("vertex_indexer").get<std::string>()
        );

        loaded.cell_state_table =
            state_table_from_json(j.at("cell_state_table"));

        loaded.vertex_table =
            state_table_row_from_json(j.at("vertex_table"));

        return loaded;

    } catch (const std::exception& ex) {
        return {};
    }

}