#!/usr/bin/env bash

set -euo pipefail

readonly TTK_REVISION="f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0"
readonly REPOSITORY_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
work_root_input="${1:-${REPOSITORY_ROOT}/../work/ttk-benchmark}"
mkdir -p "${work_root_input}"
readonly WORK_ROOT="$(cd "${work_root_input}" && pwd)"
readonly TTK_SOURCE="${WORK_ROOT}/ttk-${TTK_REVISION}"
readonly BUILD_DIR="${WORK_ROOT}/build-${TTK_REVISION}"
readonly BENCHMARK_TARGET="${2:-morseframes_ttk_gradient_benchmark}"
case "${BENCHMARK_TARGET}" in
  morseframes_ttk_gradient_benchmark|morseframes_resident_gradient_benchmark) ;;
  *) printf 'Unknown benchmark target: %s\n' "${BENCHMARK_TARGET}" >&2; exit 1 ;;
esac

if [[ ! -d "${TTK_SOURCE}/.git" ]]; then
  git init "${TTK_SOURCE}" >&2
  git -C "${TTK_SOURCE}" remote add origin \
    https://github.com/topology-tool-kit/ttk.git
fi
git -C "${TTK_SOURCE}" fetch --depth 1 origin "${TTK_REVISION}" >&2
git -C "${TTK_SOURCE}" checkout --detach "${TTK_REVISION}" >&2

cmake_args=(
  -S "${REPOSITORY_ROOT}/benchmarks/ttk"
  -B "${BUILD_DIR}"
  -DCMAKE_BUILD_TYPE=Release
  -DTTK_SOURCE_DIR="${TTK_SOURCE}"
  -DTTK_ENABLE_OPENMP=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_ZFP=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_Torch=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_CGAL=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_EMBREE=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_Graphviz=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_TBB=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON
)

if [[ "$(uname -s)" == "Darwin" ]]; then
  mac_arch="$(uname -m)"
  if [[ "$(sysctl -n hw.optional.arm64 2>/dev/null || true)" == "1" ]]; then
    mac_arch="arm64"
  fi
  cmake_args+=("-DCMAKE_OSX_ARCHITECTURES=${mac_arch}")
  libomp_root=""
  for candidate in /opt/homebrew/opt/libomp /usr/local/opt/libomp; do
    if [[ -f "${candidate}/lib/libomp.dylib" ]] &&
       file "${candidate}/lib/libomp.dylib" | grep -q "${mac_arch}"; then
      libomp_root="${candidate}"
      break
    fi
  done
  if [[ -n "${libomp_root}" ]]; then
    cmake_args+=(
      "-DOpenMP_C_FLAGS=-Xpreprocessor -fopenmp"
      "-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp"
      -DOpenMP_C_LIB_NAMES=omp
      -DOpenMP_CXX_LIB_NAMES=omp
      "-DOpenMP_omp_LIBRARY=${libomp_root}/lib/libomp.dylib"
      "-DOpenMP_C_INCLUDE_DIR=${libomp_root}/include"
      "-DOpenMP_CXX_INCLUDE_DIR=${libomp_root}/include"
    )
  fi
fi

cmake "${cmake_args[@]}" >&2
cmake --build "${BUILD_DIR}" --target "${BENCHMARK_TARGET}" --parallel >&2
printf '%s\n' "${BUILD_DIR}/${BENCHMARK_TARGET}"
