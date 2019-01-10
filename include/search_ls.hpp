#pragma once

#include "Park.hpp"

#include "nrproblem.hpp"
#include "error_norm.hpp"
#include "options.hpp"

namespace PaNr {

/* main search function */
template <typename GT,
          typename LT,
          typename FT,
          typename ST,
          typename DST,
          typename X,
          std::size_t NBORDER>
void _search_ls(GT&&                 G, // nonlinear propagator
                LT&&                 L, // linearised propagator
                ST&&                 S, // tuple of shift operators (can be empty tuple)
                DT&&                 D, //
                DVector<X, NBORDER>& z,
                Options              opt) {

    // define object for solving the NR iterations, coupling it
    // to z. Every change in z is seen by nrp.
    auto nrp = NRProblem(G, L, S, D, z);

    // temporary object for some calculations
    auto tmp = z;

    // init  parameters
    int iter = 0;

    // calc initial norm
    auto e_norm = error_norm(G, S, z, tmp);

    // main looping
    while (true) {
        // update the linear problem with the current solution
        nrp.update();

        // solve linear problem in place, overwriting the rhs
        solve_gmres(nrp.lhs(),
                    nrp.rhs(),
                    opts.gmres_rtol,
                    opts.gmres_maxiter,
                    opts.gmres_verbose);

        // perform line search with backtracking
        lambda, e_norm = linesearch(G, S, z, nrp.rhs());

        // add correction
        z = z + lambda * rnp.rhs();

        // print output

        // check tolerances
        if (iter >= opts.maxiter ||              //
            res_norm <= opts.res_norm_tol ||     //
            norm(nrp.rhs()) <= opts.dz_norm_tol) //
            break;

        // still work to do :(
        iter++;
    }
}
}