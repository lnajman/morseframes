/*    Candidate example for the Gudhi Library - https://gudhi.inria.fr/
 *    The Gudhi Library is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Laurent Najman
 *
 *    Copyright (C) 2026
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Morse_persistence.h>
#include <gudhi/Simplex_tree.h>

#include <algorithm>
#include <iostream>
#include <vector>

using Simplex_tree = Gudhi::Simplex_tree<>;
using Morse_sequence_strategy = Gudhi::morse_persistence::Morse_sequence_strategy;

void insert_simplex(Simplex_tree& simplex_tree,
                    std::initializer_list<Simplex_tree::Vertex_handle> vertices,
                    Simplex_tree::Filtration_value filtration) {
  std::vector<Simplex_tree::Vertex_handle> simplex(vertices);
  std::sort(simplex.begin(), simplex.end());
  simplex_tree.insert_simplex_and_subfaces(simplex, filtration);
}

void print_simplex(const Simplex_tree& simplex_tree, Simplex_tree::Simplex_handle simplex) {
  std::cout << '{';
  bool first = true;
  for (auto vertex : simplex_tree.simplex_vertex_range(simplex)) {
    if (!first) {
      std::cout << ',';
    }
    std::cout << vertex;
    first = false;
  }
  std::cout << '}';
}

int main() {
  Simplex_tree simplex_tree;

  // A filled triangle is inserted on one plateau, then a tail edge appears later.
  // Morse persistence handles the plateau directly, without lower-star refinement.
  insert_simplex(simplex_tree, {0}, 0.);
  insert_simplex(simplex_tree, {1}, 0.);
  insert_simplex(simplex_tree, {2}, 0.);
  insert_simplex(simplex_tree, {3}, 0.);
  insert_simplex(simplex_tree, {0, 1}, 0.);
  insert_simplex(simplex_tree, {0, 2}, 0.);
  insert_simplex(simplex_tree, {1, 2}, 0.);
  insert_simplex(simplex_tree, {0, 1, 2}, 0.);
  insert_simplex(simplex_tree, {2, 3}, 1.);
  simplex_tree.initialize_filtration();

  auto result = Gudhi::morse_persistence::compute_morse_persistence(
      simplex_tree,
      Morse_sequence_strategy::F_MAX);

  std::cout << "Number of simplices: " << simplex_tree.num_simplices() << '\n';
  std::cout << "Number of critical simplices: " << result.sequence.critical_simplices().size() << '\n';

  std::cout << "Finite off-diagonal intervals:\n";
  for (const auto& interval : result.off_diagonal_intervals()) {
    std::cout << "  H" << interval.dimension << " [" << interval.birth_value
              << ", " << interval.death_value << ")  birth simplex ";
    print_simplex(simplex_tree, result.simplex_tree_handle(interval.birth));
    std::cout << "  death simplex ";
    print_simplex(simplex_tree, result.simplex_tree_handle(interval.death));
    std::cout << '\n';
  }

  std::cout << "Essential intervals:\n";
  for (const auto& interval : result.essential_intervals()) {
    std::cout << "  H" << interval.dimension << " [" << interval.birth_value
              << ", inf)  birth simplex ";
    print_simplex(simplex_tree, result.simplex_tree_handle(interval.birth));
    std::cout << '\n';
  }

  return 0;
}

