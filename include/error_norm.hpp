#pragma once

#include "Park.hpp"

#include "nrproblem.hpp"
#include "options.hpp"

namespace PaNr {

template <typename GT,
          typename ST,
          typename X,
          std::size_t NBORDER>
double error_norm_lambda(GT&                  G,
                         ST&                  S,
                         double               lambda,
                         DVector<X, NBORDER>& z,
                         DVector<X, NBORDER>& y,
                         DVector<X, NBORDER>& tmp) {
    // understand where we are
    auto i = z.dinfo().this_rank();
    auto N = z.dinfo().comm_size();

    // Set initial state. This also performs the arithmetic on the tail
    tmp = z + lambda * y;

    // propagate by T/N
    G(tmp.head(), 0, tmp.tail(0) / N);

    // we need to shift the last segment (if we have a shift)
    if (i == N - 1 && NBORDER > 1)
        S(tmp.head(), tmp.tail(1));

    // now set tail to zero, to avoid spoling the error calculations
    err.tail() = 0;

    return norm(err);
}

template <typename GT,
          typename ST,
          typename X,
          std::size_t NBORDER>
double error_norm(GT&                  G,
                  ST&                  S,
                  DVector<X, NBORDER>& z,
                  DVector<X, NBORDER>& tmp) {
    return error_norm_lambda(G, S, 0.0, z, z, tmp);
}

}