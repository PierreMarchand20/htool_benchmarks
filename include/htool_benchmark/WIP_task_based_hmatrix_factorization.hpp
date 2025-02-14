#ifndef HTOOL_WIP_TASK_BASED_HMATRIX_FACTORIZATION_HPP
#define HTOOL_WIP_TASK_BASED_HMATRIX_FACTORIZATION_HPP

namespace htool {

struct TreeCounts {
    int call_count = 0;
    // int node_count = 0;
};

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::size_t cost_function(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {
    std::size_t nb_rows = hmatrix.get_target_cluster().get_size();
    std::size_t nb_cols = hmatrix.get_source_cluster().get_size();
    return std::size_t(nb_rows * nb_cols);
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> count_nodes(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const double criterion) { // count the number of nodes in one layer depending on criterion
    std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> result;

    if (cost_function(hmatrix) <= criterion || hmatrix.is_leaf()) {
        result.push_back(&hmatrix);
    } else {
        for (const auto &child : hmatrix.get_children()) { // postorder tree traversal
            auto local_result = count_nodes(*child.get(), criterion);
            result.insert(result.end(), local_result.begin(), local_result.end());
        }
    }
    return result;
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> find_l0(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const size_t nb_nodes_max) {
    // initialize criterion
    double criterion = cost_function(hmatrix);
    std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> old_result, result = count_nodes(hmatrix, criterion);

    if (result.size() > nb_nodes_max) {
        std::cerr << "Error: no L0 can be defined." << std::endl;
        return std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *>();

    } else {
        // dychotomy search for optimal criterion

        do {
            // save old iteration
            old_result = result;

            // if nb_nodes_max can't be reached nor surpassed we need to break loop
            if (all_leaves(result)) {
                return result;
            }
            // update criterion
            criterion /= 2;
            result = count_nodes(hmatrix, criterion);
        } while (result.size() <= nb_nodes_max);
        // {
        //     // save old iteration
        //     old_criterion = criterion;

        //     // if nb_nodes_max can't be reached nor surpassed we need to break loop
        //     if (all_leaves(count_nodes(hmatrix, criterion))) {
        //         return count_nodes(hmatrix, criterion);
        //     }

        // }
    }
    return old_result;
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
bool all_leaves(const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0) {
    return std::all_of(L0.begin(), L0.end(), [](const HMatrix<CoefficientPrecision, CoordinatePrecision> *hmatrix) {
        return hmatrix->is_leaf();
    });
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
void view_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {

    // Create a DOT file
    std::ofstream dotFile("block_tree.dot");

    // Write the DOT file contents
    dotFile << "digraph {\n";
    dotFile << "  tooltip=\"block tree\";\n";

    // Call the function to create the tree
    // create_block_tree(hmatrix, dotFile);
    TreeCounts counts;
    create_block_tree(hmatrix, dotFile, counts);

    // End of the DOT file contents
    dotFile << "}\n";

    // Close the DOT file
    dotFile.close();

    // Execute the dot command to generate an image file
    int result = system("dot -Tsvg block_tree.dot -o block_tree.svg"); // int result = in order to avoid warning
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::string get_hmatrix_id(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {
    std::string target_min = std::to_string(hmatrix.get_target_cluster().get_offset());
    std::string target_max = std::to_string(hmatrix.get_target_cluster().get_offset() + hmatrix.get_target_cluster().get_size() - 1);
    std::string source_min = std::to_string(hmatrix.get_source_cluster().get_offset());
    std::string source_max = std::to_string(hmatrix.get_source_cluster().get_offset() + hmatrix.get_source_cluster().get_size() - 1);

    return target_min + "_" + target_max + "_" + source_min + "_" + source_max;
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
void create_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, std::ostream &dotFile, TreeCounts &counts) {
    counts.call_count++;

    // Add root
    if (counts.call_count == 1) {
        auto hmatrix_info = get_hmatrix_information(hmatrix);
        dotFile << "    H_" << get_hmatrix_id(hmatrix) << " [tooltip=\"";
        for (const auto &info : hmatrix_info) {
            dotFile << info.first << ": " << info.second << "\\n";
        }
        dotFile << "\"];\n";
    }

    // Add child nodes
    for (const auto &child : hmatrix.get_children()) {
        auto child_info = get_hmatrix_information(*child.get());
        dotFile << "    H_" << get_hmatrix_id(*child.get()) << " [tooltip=\"";
        for (const auto &info : child_info) {
            dotFile << info.first << ": " << info.second << "\\n";
        }
        dotFile << "\"];\n";
        dotFile << "    H_" << get_hmatrix_id(hmatrix) << " -> H_" << get_hmatrix_id(*child.get()) << ";\n";

        // Recursively add child nodes
        create_block_tree(*child.get(), dotFile, counts);
    }
}

// template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
// void create_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, std::ostream &dotFile, TreeCounts &counts) {
//     counts.call_count++;

//     // Add root
//     if (counts.call_count == 1) {
//         auto hmatrix_info = get_hmatrix_information(hmatrix);
//         dotFile << "    H0" << " [tooltip=\"";
//         for (const auto &info : hmatrix_info) {
//             dotFile << info.first << ": " << info.second << "\\n";
//         }
//         dotFile << "\"];\n";
//         counts.node_count++;
//     }

//     // Add child nodes
//     int local_node_count = 0;
//     for (const auto &child : hmatrix.get_children()) {
//         dotFile << "    H" << counts.node_count << " [tooltip=\"";
//         auto child_info = get_hmatrix_information(*child.get());
//         for (const auto &info : child_info) {
//             dotFile << info.first << ": " << info.second << "\\n";
//         }
//         local_node_count++;
//         counts.node_count++;
//         dotFile << "local_node_count = " << local_node_count << "\\n";
//         dotFile << "node_count = " << counts.node_count << "\\n";
//         dotFile << "\"];\n";
//         dotFile << "    H" << counts.call_count - local_node_count << " -> H" << counts.call_count << ";\n";

//         // Recursively add child nodes
//         create_block_tree(*child.get(), dotFile, counts);
//         // counts.call_count++;
//     }
// }

} // namespace htool

#endif