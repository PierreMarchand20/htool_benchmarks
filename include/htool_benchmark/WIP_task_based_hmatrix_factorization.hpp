#ifndef HTOOL_WIP_TASK_BASED_HMATRIX_FACTORIZATION_HPP
#define HTOOL_WIP_TASK_BASED_HMATRIX_FACTORIZATION_HPP
/*
This file contain usefull functions for task based hmatrix factorization as described in the phD thesis of Benoit Lizé {https://theses.hal.science/tel-01244260}.
*/

namespace htool {

struct TreeCounts {
    size_t call_count = 0;
    size_t node_count = 0;
};

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::size_t cost_function(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {
    /*
    The cost_function associates with a node of the group tree a score representing an estimate of the amount of work associated with this leaf. Below is the naive expression of this score.
    */
    std::size_t nb_rows = hmatrix.get_target_cluster().get_size();
    std::size_t nb_cols = hmatrix.get_source_cluster().get_size();
    return std::size_t(nb_rows * nb_cols);
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> count_nodes(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const double criterion) {
    /*
    This function performs a postorder tree traversal of the group tree and returns a vector of all the nodes of the group tree that have a cost less than criterion.

    The cost of a node is given by the cost_function: it is the product of the number of points in the target cluster and the number of points in the source cluster.

    The criterion is the maximum cost of the nodes to be returned.

    The result is a vector of pointers to the nodes of the group tree that have a cost less than criterion.
    */
    std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> result;

    if (cost_function(hmatrix) <= criterion || hmatrix.is_leaf()) {
        // if the node is a leaf or its cost is less than criterion, add it to the result
        result.push_back(&hmatrix);
    } else {
        // if the node is not a leaf, traverse its children
        for (const auto &child : hmatrix.get_children()) {
            // perform a postorder tree traversal of the subtree rooted at child
            auto local_result = count_nodes(*child.get(), criterion);
            // add the result of the subtree traversal to the result
            result.insert(result.end(), local_result.begin(), local_result.end());
        }
    }
    return result;
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> find_l0(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const size_t nb_nodes_max) {
    /*
    This function returns a vector of all the nodes in the group tree of the HMatrix
    that have a cost less than an optimal criterion based on the maximum number of nodes
    allowed (nb_nodes_max).

    The criterion starts as the cost of the root node and is adjusted through a dichotomy
    search to find an optimal set of nodes that does not exceed nb_nodes_max.
    */

    // Initialize criterion with the cost of the root node
    double criterion = cost_function(hmatrix);

    // Find initial nodes that meet the criterion
    std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> old_result, result = count_nodes(hmatrix, criterion);

    // Check if the initial result exceeds the maximum allowed nodes
    if (result.size() > nb_nodes_max) {
        std::cerr << "Error: no L0 can be defined." << std::endl;
        return std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *>();
    } else {
        // Perform a dichotomy search to find the optimal criterion
        do {
            // Save the result of this iteration
            old_result = result;

            // If all nodes are leaves, return the result as it cannot be further divided
            if (all_leaves(result)) {
                return result;
            }

            // Reduce the criterion to explore smaller cost nodes
            criterion /= 2;

            // Update the result with the new criterion
            result = count_nodes(hmatrix, criterion);
        } while (result.size() <= nb_nodes_max); // Loop until the result size exceeds nb_nodes_max
    }

    // Return the last valid result
    return old_result;
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
bool all_leaves(const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0) {
    // This function checks if all elements in the vector L0 are leaf nodes.

    // Use std::all_of to determine if all elements satisfy the condition
    return std::all_of(L0.begin(), L0.end(), [](const HMatrix<CoefficientPrecision, CoordinatePrecision> *hmatrix) {
        // Check if the current HMatrix is a leaf
        return hmatrix->is_leaf();
    });
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
void view_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {
    /*
    This function visualizes the block tree of an HMatrix by generating a DOT file and converting it into an SVG image.
    The DOT file captures the hierarchical structure of the HMatrix.
    */
    // Find the initial set of nodes L0 with a maximum size of 64
    std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> L0 = find_l0(hmatrix, 64);

    // Create and open a DOT file to write the block tree data
    std::ofstream dotFile("block_tree.dot");

    // Start the DOT file content
    dotFile << "digraph {\n";

    // Create the block tree by adding nodes and edges to the DOT file
    TreeCounts counts;
    create_block_tree(hmatrix, dotFile, counts, L0);

    // Add tooltip information to the DOT file summarizing the block tree
    dotFile << "  tooltip=\"Block tree information: \\n"
            << "number of nodes: " << counts.node_count << "\\n"
            << "number of nodes in L0: " << L0.size() << "\";\n";
    dotFile << "}\n";

    // Close the DOT file after writing all content
    dotFile.close();

    // Execute the DOT command to convert the DOT file into an SVG image
    // The result is stored in an SVG file named 'block_tree.svg'
    int result = system("dot -Tsvg block_tree.dot -o block_tree.svg"); // int result = to avoid warning about unused return value
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::string get_hmatrix_id(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {
    // Get the minimum and maximum offsets of the target and source clusters
    std::string target_min = std::to_string(hmatrix.get_target_cluster().get_offset());
    std::string target_max = std::to_string(hmatrix.get_target_cluster().get_offset() + hmatrix.get_target_cluster().get_size() - 1);
    std::string source_min = std::to_string(hmatrix.get_source_cluster().get_offset());
    std::string source_max = std::to_string(hmatrix.get_source_cluster().get_offset() + hmatrix.get_source_cluster().get_size() - 1);

    // Create the identifier by concatenating the offsets
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
    /*
    This function adds a node to the block tree.
    It creates a new node with the given hmatrix and adds it to the dot file.
    It also adds the node information in a tooltip.
    If the node is in L0, it colors it in light blue.
    The function also increments the node count.
    */

    // Define node
    dotFile << "    H_" << get_hmatrix_id(hmatrix) << " [tooltip=\"Node information: \\n";

    // Add the node information in a tooltip.
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
    // Increment the call count
    counts.call_count++;

    // Add root node
    if (counts.call_count == 1) {
        // Add the root node to the dot file
        add_node_to_block_tree(hmatrix, dotFile, counts, L0);
    }

    // Add child nodes
    for (const auto &child : hmatrix.get_children()) {
        // Add the child node to the dot file
        add_node_to_block_tree(*child.get(), dotFile, counts, L0);

        // Add an edge between the parent and child nodes
        dotFile << "    H_" << get_hmatrix_id(hmatrix) << " -> H_" << get_hmatrix_id(*child.get()) << " [tooltip=\"" << get_hmatrix_id(hmatrix) << " -> " << get_hmatrix_id(*child.get()) << "\"];\n";

        // Recursively add child nodes
        create_block_tree(*child.get(), dotFile, counts, L0);
    }
}

template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> enumerate_dependances(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0, const HMatrix<CoefficientPrecision, CoordinatePrecision> &root) {
    /*
    Let be M a node of a block tree. We want to enumerate all the nodes that M is dependant on and are in L0. These nodes are called dependances of M and are noted D(M). D(M) is defined in function of the position of M vis a vis L0. M can be on, above or below L0.
    - if M is on L0, then D(M) = {M}
    - if M is above L0, then D(M) = {L0 & child_of(M)}
    - if M is below L0, then D(M) = {parent_of(M) & L0}
    */

    // Case 1 : hmatrix is in L0
    if (std::find(L0.begin(), L0.end(), &hmatrix) != L0.end()) {
        return {&hmatrix};
    }

    // Case 2 : hmatrix is above L0. Find descendants of hmatrix in L0
    std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> result;
    for (const auto &child : hmatrix.get_children()) {
        if (std::find(L0.begin(), L0.end(), child.get()) != L0.end()) {
            result.push_back(child.get());
        }
        auto dependances = enumerate_dependances(*child.get(), L0);
        result.insert(result.end(), dependances.begin(), dependances.end());
    }
    return result;

    // Case 3 : hmatrix is below L0. Find the ancestor of hmatrix in L0
    return get_parent(hmatrix, root);
}

template <typename CoefficientPrecision, typename CoordinatePrecision>
const HMatrix<CoefficientPrecision, CoordinatePrecision> *get_parent(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const HMatrix<CoefficientPrecision, CoordinatePrecision> &root) {
    if (root.get_children().empty()) {
        return nullptr;
    }
    for (const auto &child : root.get_children()) {
        if (child.get() == &hmatrix) {
            return &root;
        }
        const HMatrix<CoefficientPrecision, CoordinatePrecision> *parent = get_parent(hmatrix, *child.get());
        if (parent != nullptr) {
            return parent;
        }
    }
    return nullptr;
}

} // namespace htool

#endif