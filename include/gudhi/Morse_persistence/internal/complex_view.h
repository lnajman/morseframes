#pragma once

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

#include "gudhi/Morse_persistence/internal/filtered_complex.h"

namespace Gudhi { namespace morse_persistence { namespace internal {

namespace detail {

template <class View, class = void>
struct has_complex_view_api : std::false_type {};

template <class Range, class Id, class = void>
struct is_id_range : std::false_type {};

template <class Range, class Id>
struct is_id_range<
    Range,
    Id,
    std::void_t<
        decltype(std::declval<const std::remove_reference_t<Range>&>().size()),
        decltype(std::begin(std::declval<const std::remove_reference_t<Range>&>())),
        decltype(std::end(std::declval<const std::remove_reference_t<Range>&>())),
        decltype(*std::begin(std::declval<const std::remove_reference_t<Range>&>()))>>
    : std::conjunction<
          std::is_convertible<
              decltype(std::declval<const std::remove_reference_t<Range>&>().size()),
              std::size_t>,
          std::is_convertible<
              decltype(*std::begin(std::declval<const std::remove_reference_t<Range>&>())),
              Id>> {};

template <class View>
struct has_complex_view_api<
    View,
    std::void_t<
        decltype(std::declval<const View&>().size()),
        decltype(std::declval<const View&>().dimension(std::declval<SimplexId>())),
        decltype(std::declval<const View&>().level(std::declval<SimplexId>())),
        decltype(std::declval<const View&>().filtration(std::declval<SimplexId>())),
        decltype(std::declval<const View&>().vertices(std::declval<SimplexId>())),
        decltype(std::declval<const View&>().boundary(std::declval<SimplexId>())),
        decltype(std::declval<const View&>().coboundary(std::declval<SimplexId>())),
        decltype(std::declval<const View&>().filtration_order()),
        decltype(std::declval<const View&>().simplices_of_level(std::declval<LevelId>())),
        decltype(std::declval<const View&>().level_values()),
        decltype(std::declval<const View&>().num_levels())>>
    : std::conjunction<
          std::is_convertible<decltype(std::declval<const View&>().size()), std::size_t>,
          std::is_convertible<
              decltype(std::declval<const View&>().dimension(std::declval<SimplexId>())),
              std::uint16_t>,
          std::is_convertible<
              decltype(std::declval<const View&>().level(std::declval<SimplexId>())),
              LevelId>,
          std::is_convertible<
              decltype(std::declval<const View&>().filtration(std::declval<SimplexId>())),
              double>,
          is_id_range<
              decltype(std::declval<const View&>().vertices(std::declval<SimplexId>())),
              VertexId>,
          is_id_range<
              decltype(std::declval<const View&>().boundary(std::declval<SimplexId>())),
              SimplexId>,
          is_id_range<
              decltype(std::declval<const View&>().coboundary(std::declval<SimplexId>())),
              SimplexId>,
          std::is_convertible<decltype(std::declval<const View&>().filtration_order()),
                              const std::vector<SimplexId>&>,
          std::is_convertible<
              decltype(std::declval<const View&>().simplices_of_level(std::declval<LevelId>())),
              const std::vector<SimplexId>&>,
          std::is_convertible<decltype(std::declval<const View&>().level_values()),
                              const std::vector<double>&>,
          std::is_convertible<decltype(std::declval<const View&>().num_levels()),
                              std::size_t>> {};

}  // namespace detail

template <class View>
struct is_complex_view : detail::has_complex_view_api<View> {};

template <class View>
constexpr bool is_complex_view_v = is_complex_view<View>::value;

template <class View>
void require_complex_view() {
  static_assert(
      is_complex_view_v<View>,
      "A Morse complex view must expose size(), dimension(id), level(id), "
      "filtration(id), vertices(id), boundary(id), coboundary(id), "
      "filtration_order(), simplices_of_level(level), level_values(), and num_levels().");
}

}  // namespace internal
}  // namespace morse_persistence
}  // namespace Gudhi
