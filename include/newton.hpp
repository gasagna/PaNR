#pragma once

#include <algorithm>
#include <tuple>

#include "ParK"

#include "nrproblem.hpp"
#include "options.hpp"

namespace PaNr {

template <typename GT,
          typename LT,
          typename FT,
          typename X>
void search(GT&            G, // nonlinear state-transition propagator
            LT&            L, // linearised state-transition propagator
            FT&            F, // right hand side of governing equations
            DVector<X, 1>& z, // initial guess with no shifts
            Options        opts = Options()) {
    
    opts.method == "linesearch"
        ? _search_ls(G, L, std::tie(std::nullptr), std::tie(F), z, opts)
        : _search_hs(G, L, std::tie(std::nullptr), std::tie(F), z, opts);
}

template <typename GT,
          typename LT,
          typename FT,
          typename ST,
          typename DST,
          typename X>
void search(GT&            G, // nonlinear propagator
            LT&            L, // linearised propagator
            FT&            F, // right hand side
            ST&            S, // tuple of shift operators (can be empty tuple)
            DST&           D, // tuple of derivatives of the shift operators
            DVector<X, 2>& z,
            Options        opts = Options()) {
    opts.method == "linesearch"
        ? _search_ls(G, L, S, std::tie(F, DS), z, opts)
        : _search_hs(G, L, S, std::tie(F, DS), z, opts);
}

}