#ifndef HTOOL_WIP_TASK_BASED_HMATRIX_FACTORIZATION_HPP
#define HTOOL_WIP_TASK_BASED_HMATRIX_FACTORIZATION_HPP

namespace htool {

struct TreeCounts {
    size_t call_count = 0;
    size_t node_count = 0;
};

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::size_t cost_function(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {
    std::size_t nb_rows = hmatrix.get_target_cluster().get_size();
    std::size_t nb_cols = hmatrix.get_source_cluster().get_size();
    return std::size_t(nb_rows * nb_cols);
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> count_nodes(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const double criterion) {

    // count the number of nodes in one layer depending on criterion
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
        // dichotomy search for optimal criterion
        do {
            // save old iteration
            old_result = result;

            // if nb_nodes_max can't be reached nor surpassed we need to break loop
            if (all_leaves(result)) {
                return result;
            }

            // update criterion
            criterion /= 2;

            // update result
            result = count_nodes(hmatrix, criterion);
        } while (result.size() <= nb_nodes_max);
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
    // Find L0
    std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> L0 = find_l0(hmatrix, 64);

    // Create a DOT file
    std::ofstream dotFile("block_tree.dot");

    // Write the DOT file contents
    dotFile << "digraph {\n";

    // Call the function to create the tree
    TreeCounts counts;
    create_block_tree(hmatrix, dotFile, counts, L0);

    // End of the DOT file contents
    dotFile << "  tooltip=\"Block tree information: \\n"
            << "number of nodes: " << counts.node_count << "\\n"
            << "number of nodes in L0: " << L0.size() << "\";\n";
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

// template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
// void create_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, std::ostream &dotFile, TreeCounts &counts, const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0) {
//     counts.call_count++;

//     // Add root
//     if (counts.call_count == 1) {
//         // Define node
//         dotFile << "    H_" << get_hmatrix_id(hmatrix) << " [tooltip=\"Node information: \\n";

//         // Add node information
//         auto hmatrix_info = get_hmatrix_information(hmatrix);
//         for (const auto &info : hmatrix_info) {
//             dotFile << info.first << ": " << info.second << "\\n";
//         }

//         // Check if the current node is in L0 and color it in light blue if true
//         dotFile << "\"";
//         if (std::find(L0.begin(), L0.end(), &hmatrix) != L0.end()) {
//             dotFile << ", style=filled, fillcolor=lightblue";
//         }

//         // End of node definition
//         dotFile << "];\n";

//         // Increment node count
//         counts.node_count++;
//     }

//     // Add child nodes
//     for (const auto &child : hmatrix.get_children()) {
//         // Define node
//         dotFile << "    H_" << get_hmatrix_id(*child.get()) << " [tooltip=\"Node information: \\n";

//         // Add node information
//         auto child_info = get_hmatrix_information(*child.get());
//         for (const auto &info : child_info) {
//             dotFile << info.first << ": " << info.second << "\\n";
//         }

//         // Check if the child node is in L0 and color it in light blue if true
//         dotFile << "\"";
//         if (std::find(L0.begin(), L0.end(), child.get()) != L0.end()) {
//             dotFile << ", style=filled, fillcolor=lightblue";
//         }

//         // End of node definition
//         dotFile << "];\n";

//         // Add edge
//         dotFile << "    H_" << get_hmatrix_id(hmatrix) << " -> H_" << get_hmatrix_id(*child.get()) << " [tooltip=\"" << get_hmatrix_id(hmatrix) << " -> " << get_hmatrix_id(*child.get()) << "\"];\n";

//         // Increment node count
//         counts.node_count++;

//         // Recursively add child nodes
//         create_block_tree(*child.get(), dotFile, counts, L0);
//     }
// }

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
void add_node_to_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, std::ostream &dotFile, TreeCounts &counts, const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0) {
    // Define node
    dotFile << "    H_" << get_hmatrix_id(hmatrix) << " [tooltip=\"Node information: \\n";

    // Add node information
    auto hmatrix_info = get_hmatrix_information(hmatrix);
    for (const auto &info : hmatrix_info) {
        dotFile << info.first << ": " << info.second << "\\n";
    }

    // Check if the current node is in L0 and color it in light blue if true
    dotFile << "\"";
    if (std::find(L0.begin(), L0.end(), &hmatrix) != L0.end()) {
        dotFile << ", style=filled, fillcolor=lightblue";
    }

    // End of node definition
    dotFile << "];\n";

    // Increment node count
    counts.node_count++;
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
void create_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, std::ostream &dotFile, TreeCounts &counts, const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0) {
    counts.call_count++;

    // Add root
    if (counts.call_count == 1) {
        add_node_to_block_tree(hmatrix, dotFile, counts, L0);
    }

    // Add child nodes
    for (const auto &child : hmatrix.get_children()) {
        add_node_to_block_tree(*child.get(), dotFile, counts, L0);

        // Add edge
        dotFile << "    H_" << get_hmatrix_id(hmatrix) << " -> H_" << get_hmatrix_id(*child.get()) << " [tooltip=\"" << get_hmatrix_id(hmatrix) << " -> " << get_hmatrix_id(*child.get()) << "\"];\n";

        // Recursively add child nodes
        create_block_tree(*child.get(), dotFile, counts, L0);
    }
}

} // namespace htool

#endif