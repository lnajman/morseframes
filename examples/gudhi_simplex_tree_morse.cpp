#include <algorithm>
#include <iostream>
#include <vector>

#include <gudhi/Morse_persistence.h>
#include <gudhi/Simplex_tree.h>

namespace {

using SimplexTree = Gudhi::Simplex_tree<>;
namespace mp = Gudhi::morse_persistence;

void insert(SimplexTree& simplex_tree,
            std::initializer_list<SimplexTree::Vertex_handle> vertices,
            double filtration) {
  std::vector<SimplexTree::Vertex_handle> simplex(vertices);
  std::sort(simplex.begin(), simplex.end());
  simplex_tree.insert_simplex_and_subfaces(simplex, filtration);
}

SimplexTree make_example_simplex_tree() {
  SimplexTree simplex_tree;

  insert(simplex_tree, {0}, 0.0);
  insert(simplex_tree, {1}, 0.0);
  insert(simplex_tree, {2}, 0.0);
  insert(simplex_tree, {3}, 0.0);
  insert(simplex_tree, {0, 1}, 0.0);
  insert(simplex_tree, {0, 2}, 0.0);
  insert(simplex_tree, {1, 2}, 0.0);
  insert(simplex_tree, {0, 1, 2}, 0.0);
  insert(simplex_tree, {2, 3}, 1.0);

  simplex_tree.initialize_filtration();
  return simplex_tree;
}

void print_simplex(const SimplexTree& simplex_tree, SimplexTree::Simplex_handle handle) {
  std::cout << '{';
  bool first = true;
  for (auto vertex : simplex_tree.simplex_vertex_range(handle)) {
    if (!first) {
      std::cout << ',';
    }
    std::cout << vertex;
    first = false;
  }
  std::cout << '}';
}

}  // namespace

int main() {
  auto simplex_tree = make_example_simplex_tree();
  constexpr auto strategy = mp::Morse_sequence_strategy::F_MAX;
  const auto result = mp::compute_morse_persistence(
      simplex_tree,
      strategy);

  std::cout << "strategy\t" << mp::strategy_name(strategy) << '\n';
  std::cout << "simplices\t" << result.view.size() << '\n';
  std::cout << "critical_simplices\t" << result.sequence.critical_simplices().size()
            << '\n';

  for (const auto& pair : result.off_diagonal_intervals()) {
    std::cout << "finite\t" << pair.dimension << '\t' << pair.birth_value << '\t'
              << pair.death_value << "\tbirth=";
    print_simplex(simplex_tree, result.simplex_tree_handle(pair.birth));
    std::cout << "\tdeath=";
    print_simplex(simplex_tree, result.simplex_tree_handle(pair.death));
    std::cout << '\n';
  }

  for (const auto& interval : result.essential_intervals()) {
    std::cout << "essential\t" << interval.dimension << '\t' << interval.birth_value
              << "\tinf\tbirth=";
    print_simplex(simplex_tree, result.simplex_tree_handle(interval.birth));
    std::cout << '\n';
  }

  return 0;
}
