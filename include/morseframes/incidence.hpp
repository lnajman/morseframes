#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "morseframes/complex_view.hpp"
#include "morseframes/field_arithmetic.hpp"

namespace morseframes {

namespace detail {

template <class View, class = void>
struct has_boundary_coefficient_api : std::false_type {};

template <class View>
struct has_boundary_coefficient_api<
    View,
    std::void_t<decltype(std::declval<const View&>().boundary_coefficient(
        std::declval<SimplexId>(),
        std::declval<std::size_t>(),
        std::declval<std::uint32_t>()))>>
    : std::is_convertible<decltype(std::declval<const View&>().boundary_coefficient(
                              std::declval<SimplexId>(),
                              std::declval<std::size_t>(),
                              std::declval<std::uint32_t>())),
                          std::uint32_t> {};

}  // namespace detail

template <class ComplexView>
inline std::uint32_t boundary_incidence_coefficient(const ComplexView& complex,
                                                    SimplexId cell,
                                                    std::size_t boundary_index,
                                                    std::uint32_t modulus) {
  static_assert(is_complex_view_v<ComplexView>,
                "boundary_incidence_coefficient requires a Morse complex-view type.");
  if constexpr (detail::has_boundary_coefficient_api<ComplexView>::value) {
    return complex.boundary_coefficient(cell, boundary_index, modulus);
  } else {
    (void)complex;
    (void)cell;
    return boundary_coefficient(boundary_index, modulus);
  }
}

template <class ComplexView>
inline std::uint32_t boundary_incidence_coefficient_of_face(
    const ComplexView& complex,
    SimplexId coface,
    SimplexId face,
    std::uint32_t modulus) {
  static_assert(is_complex_view_v<ComplexView>,
                "boundary_incidence_coefficient_of_face requires a Morse "
                "complex-view type.");
  const auto& boundary = complex.boundary(coface);
  for (std::size_t boundary_index = 0; boundary_index < boundary.size(); ++boundary_index) {
    if (boundary[boundary_index] == face) {
      return boundary_incidence_coefficient(complex, coface, boundary_index, modulus);
    }
  }
  throw std::logic_error("Expected a codimension-one face/coface incidence.");
}

}  // namespace morseframes
