add_executable_with_targets(
  Morse_persistence_test_simplex_tree
  test_morse_persistence_simplex_tree.cpp
  TBB::tbb)

add_test(
  NAME Morse_persistence_test_simplex_tree
  COMMAND $<TARGET_FILE:Morse_persistence_test_simplex_tree>)
