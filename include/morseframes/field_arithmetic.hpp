#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace morseframes {

inline bool is_prime_field_characteristic(std::uint32_t modulus) {
  if (modulus < 2) {
    return false;
  }
  if (modulus == 2) {
    return true;
  }
  if (modulus % 2 == 0) {
    return false;
  }
  for (std::uint32_t divisor = 3; divisor <= modulus / divisor; divisor += 2) {
    if (modulus % divisor == 0) {
      return false;
    }
  }
  return true;
}

inline void validate_prime_field_characteristic(std::uint32_t modulus) {
  if (!is_prime_field_characteristic(modulus)) {
    throw std::invalid_argument(
        "Coefficient modulus must be prime; composite moduli are not supported.");
  }
}

inline std::uint32_t boundary_coefficient(std::size_t removed_index,
                                          std::uint32_t modulus) {
  return removed_index % 2 == 0 ? 1u : modulus - 1u;
}

inline std::uint32_t modp_multiply(std::uint32_t lhs,
                                   std::uint32_t rhs,
                                   std::uint32_t modulus) {
  return static_cast<std::uint32_t>(
      (static_cast<std::uint64_t>(lhs) * static_cast<std::uint64_t>(rhs)) % modulus);
}

inline std::uint32_t modp_power(std::uint32_t base,
                                std::uint32_t exponent,
                                std::uint32_t modulus) {
  std::uint32_t result = 1;
  base %= modulus;
  while (exponent > 0) {
    if ((exponent & 1u) != 0u) {
      result = modp_multiply(result, base, modulus);
    }
    exponent >>= 1u;
    if (exponent != 0) {
      base = modp_multiply(base, base, modulus);
    }
  }
  return result;
}

inline std::uint32_t modp_inverse(std::uint32_t value, std::uint32_t modulus) {
  if (value % modulus == 0) {
    throw std::logic_error("Cannot invert zero modulo the coefficient field.");
  }
  return modp_power(value, modulus - 2u, modulus);
}

}  // namespace morseframes
