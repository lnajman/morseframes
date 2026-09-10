#pragma once

#include "morseframes/morse_sequence.hpp"

namespace morseframes {

// Lightweight, validated initialization for RK-only callers. Private inheritance
// reuses the existing kernel, callback, profiling and replay implementation, but
// does not expose strategies that need FSequenceBuilder's omitted caches. There
// is no mutable lazy state; the general-purpose builder remains unchanged.
template <class ComplexView = FilteredSimplicialComplex>
class ReductionKernelSequenceBuilder : private FSequenceBuilder<ComplexView> {
  using Base = FSequenceBuilder<ComplexView>;

 public:
  explicit ReductionKernelSequenceBuilder(
      const ComplexView& complex,
      MorseSequenceBuildMetrics* sequence_metrics = nullptr,
      bool detailed_reduction_kernel_metrics = true)
      : Base(complex, sequence_metrics, detailed_reduction_kernel_metrics,
             typename Base::ReductionKernelOnlyTag{}) {}

  using Base::build_flooding_reduction_kernel;
  using Base::build_flooding_reduction_kernel_parallel;
  using Base::build_flooding_reduction_kernel_with_step_callback;
  using Base::build_flooding_reduction_kernel_parallel_with_step_callback;
  using Base::build_flooding_reduction_kernel_with_execution_options;
};

template <class ComplexView>
ReductionKernelSequenceBuilder(const ComplexView&)
    -> ReductionKernelSequenceBuilder<ComplexView>;

}  // namespace morseframes
