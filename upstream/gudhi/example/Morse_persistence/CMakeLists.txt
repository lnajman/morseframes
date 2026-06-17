add_executable_with_targets(
  Morse_persistence_example_from_simplex_tree
  example_morse_persistence_from_simplex_tree.cpp
  TBB::tbb)

add_test(
  NAME Morse_persistence_example_from_simplex_tree
  COMMAND $<TARGET_FILE:Morse_persistence_example_from_simplex_tree>)
