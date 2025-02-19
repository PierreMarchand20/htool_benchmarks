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

/**
 * @brief The cost_function associates with a node of the group tree a score
 * representing an estimate of the amount of work associated with this leaf.
 *
 * The cost of a node is given by the number of points in the target cluster
 * times the number of points in the source cluster.
 *
 * @param hmatrix The input hierarchical matrix.
 * @return The naive cost of the node, i.e. the number of points in the target cluster
 * times the number of points in the source cluster.
 */
template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::size_t cost_function(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {
    std::size_t nb_rows = hmatrix.get_target_cluster().get_size();
    std::size_t nb_cols = hmatrix.get_source_cluster().get_size();
    return std::size_t(nb_rows * nb_cols);
}

/**
 * @brief Performs a postorder tree traversal of the group tree and returns a vector of all the nodes of the group tree that have a cost less than criterion.
 *
 * The cost of a node is given by the cost_function: it is the product of the number of points in the target cluster and the number of points in the source cluster.
 *
 * The criterion is the maximum cost of the nodes to be returned.
 *
 * The result is a vector of pointers to the nodes of the group tree that have a cost less than criterion.
 *
 * @param hmatrix The input hierarchical matrix.
 * @param criterion The maximum cost of the nodes to be returned.
 * @return A vector of pointers to the nodes of the group tree that have a cost less than criterion.
 */
template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> count_nodes(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const double criterion) {
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

/**
 * @brief This function returns a vector of all the nodes in the group tree of the HMatrix
 * that have a cost less than an optimal criterion based on the maximum number of nodes
 * allowed (nb_nodes_max).
 *
 * The criterion starts as the cost of the root node and is adjusted through a dichotomy
 * search to find an optimal set of nodes that does not exceed nb_nodes_max.
 *
 * @param hmatrix The input hierarchical matrix.
 * @param nb_nodes_max The maximum number of nodes in the result.
 * @return A vector of all the nodes in the group tree of the HMatrix that have a cost less than an optimal criterion.
 */
template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> find_l0(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const size_t nb_nodes_max) {
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

/**
 * @brief Checks if all HMatrix nodes in the provided vector are leaf nodes.
 *
 * This function iterates through the vector L0 and determines if each
 * HMatrix node is a leaf node. It returns true if all nodes are leaves,
 * otherwise false.
 *
 * @param L0 A vector of pointers to HMatrix nodes to be checked.
 * @return True if all nodes in L0 are leaf nodes, false otherwise.
 */
template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
bool all_leaves(const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0) {
    // This function checks if all elements in the vector L0 are leaf nodes.

    // Use std::all_of to determine if all elements satisfy the condition
    return std::all_of(L0.begin(), L0.end(), [](const HMatrix<CoefficientPrecision, CoordinatePrecision> *hmatrix) {
        // Check if the current HMatrix is a leaf
        return hmatrix->is_leaf();
    });
}

/**
 * @brief Visualizes the block tree of a hierarchical matrix (HMatrix).
 *
 * This function generates a DOT file representing the hierarchical structure of the given HMatrix.
 * The DOT file is then converted into an SVG image to visualize the block tree.
 * The function first determines an initial set of nodes (L0) with a maximum size of 64.
 * It then writes the block tree data, including nodes and edges, to a DOT file.
 * Tooltip information summarizing the block tree, such as the number of nodes and nodes in L0,
 * is added to the DOT file. Finally, the DOT file is converted to an SVG image using the
 * system's DOT command.
 *
 * @param hmatrix The input hierarchical matrix to be visualized.
 */
template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
void view_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix) {
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

/**
 * @brief Creates a unique identifier for a given HMatrix based on its target and source cluster offsets.
 *
 * This function generates an identifier by concatenating the minimum and maximum offsets of the target and source clusters.
 * This identifier is used to label nodes in the block tree visualization of the HMatrix.
 *
 * @param hmatrix The input hierarchical matrix for which the identifier is generated.
 * @return A unique identifier for the given HMatrix.
 */
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

/**
 * @brief Adds a node to the block tree in the DOT file.
 *
 * This function adds a new node to the dot file.
 * It also adds the node information in a tooltip.
 * If the node is in L0, it colors it in light blue.
 * The function also increments the node count.
 *
 * @param hmatrix The input hierarchical matrix for which the node is created.
 * @param dotFile The output stream to write the node data to.
 * @param counts The TreeCounts object to increment the node count.
 * @param L0 The initial set of nodes with a maximum size of 64.
 */
template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
void add_node_to_block_tree(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, std::ostream &dotFile, TreeCounts &counts, const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0) {
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

/**
 * @brief Constructs a block tree visualization for a hierarchical matrix (HMatrix).
 *
 * This function generates the structure of a block tree by recursively adding nodes
 * and edges to a provided DOT file stream. It starts by adding the root node and
 * then iterates over the children of each node, adding them to the DOT file and
 * linking them with edges. The tree is constructed in a depth-first manner.
 * Each node and edge is annotated with tooltips for additional information.
 *
 * @param hmatrix The input hierarchical matrix for which the block tree is created.
 * @param dotFile The output stream to write the DOT representation of the block tree.
 * @param counts A TreeCounts object that tracks the number of calls and nodes processed.
 * @param L0 A vector of pointers to HMatrix nodes, representing an initial set of nodes.
 */
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

/**
 * @brief Enumerates the dependencies of a hierarchical matrix node within a block tree.
 *
 * This function identifies all nodes that a given matrix node (M) depends on, which are also present
 * in the initial set of nodes (L0). The dependencies (D(M)) are determined based on the position
 * of the node M with respect to L0. Specifically:
 * - If M is on L0, D(M) consists of just M itself.
 * - If M is above L0, D(M) consists of the intersection of L0 and the descendants of M.
 * - If M is below L0, D(M) consists of the intersection of L0 and the ancestors of M.
 *
 * @param hmatrix The matrix node whose dependencies are to be enumerated.
 * @param L0 A vector containing the initial set of nodes within the block tree.
 * @param root The root node of the hierarchical matrix tree.
 * @return A vector of pointers to the HMatrix nodes that are dependencies of the given node.
 */
template <typename CoefficientPrecision, typename CoordinatePrecision = underlying_type<CoefficientPrecision>>
std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> enumerate_dependances(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const std::vector<const HMatrix<CoefficientPrecision, CoordinatePrecision> *> &L0, const HMatrix<CoefficientPrecision, CoordinatePrecision> &root) {
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
    // Todo : rather than starting from root, search by starting from each node in L0
    return get_parent(hmatrix, root);
}

/**
 * @brief Finds the parent of a node in a hierarchical matrix tree.
 *
 * Recursively traverses the tree, starting from the root node, until it finds the parent of the given node.
 *
 * @param hmatrix The matrix node for which the parent is to be found.
 * @param root The root node of the hierarchical matrix tree.
 * @return A pointer to the parent node, or nullptr if the parent is not found.
 */
template <typename CoefficientPrecision, typename CoordinatePrecision>
const HMatrix<CoefficientPrecision, CoordinatePrecision> *get_parent(const HMatrix<CoefficientPrecision, CoordinatePrecision> &hmatrix, const HMatrix<CoefficientPrecision, CoordinatePrecision> &root) {
    // Case 1 : hmatrix is the root
    if (root.get_children().empty()) {
        return nullptr;
    }

    for (const auto &child : root.get_children()) {
        // Case 2 : hmatrix is a child of root
        if (child.get() == &hmatrix) {
            return &root;
        }

        // Case 3 : recursively search for the parent
        const HMatrix<CoefficientPrecision, CoordinatePrecision> *parent = get_parent(hmatrix, *child.get());

        // Case 4 : parent is found
        if (parent != nullptr) {
            return parent;
        }
    }
    return nullptr;
}

} // namespace htool

#endif