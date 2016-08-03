/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#ifndef _EXAHYPE_KERNELS_KERNEL_UTILS_H_
#define _EXAHYPE_KERNELS_KERNEL_UTILS_H_

#include "../../Peano/tarch/Assertions.h"

namespace kernels {

struct idx2 {
  idx2(int I, int J, int line = -1) : I_(I), J_(J), line_(line) {}

  int operator()(int i, int j) {
    assertion3(i < I_, i, I_, line_);
    assertion3(j < J_, j, J_, line_);
    return i * J_ + j;
  }

  const int I_, J_, line_;
};

struct idx3 {
  idx3(int I, int J, int K, int line = -1) : I_(I), J_(J), K_(K), line_(line) {}

  int operator()(int i, int j, int k) {
    assertion3(i < I_, i, I_, line_);
    assertion3(j < J_, j, J_, line_);
    assertion3(k < K_, k, K_, line_);
    return i * (J_ * K_) + j * K_ + k;
  }

  const int I_, J_, K_, line_;
};

struct idx4 {
  idx4(int I, int J, int K, int L, int line = -1)
      : I_(I), J_(J), K_(K), L_(L), line_(line) {}

  int operator()(int i, int j, int k, int l) {
    assertion3(i < I_, i, I_, line_);
    assertion3(j < J_, j, J_, line_);
    assertion3(k < K_, k, K_, line_);
    assertion3(l < L_, l, L_, line_);
    return i * (J_ * K_ * L_) + j * (K_ * L_) + k * L_ + l;
  }

  const int I_, J_, K_, L_, line_;
};

}  // namespace kernels

#endif  // _EXAHYPE_KERNELS_KERNEL_UTILS_H_
