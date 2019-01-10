#pragma once

#include "Park.hpp"

#include "error_norm.hpp"
#include "nrproblem.hpp"
#include "options.hpp"

namespace PaNr {

auto linesearch(GT&                  G,
                ST&                  S,
                DVector<X, NBORDER>& z,
                DVector<X, NBORDER>& dz,
                DVector<X, NBORDER>& tmp) {

    // start with full newton step
    double lambda = 1.0;

    // but in the first iterations it might be possible that
    // this results in a period that is too small, or even
    // negative. Hence, we limit lambda such that it will
    // only give an integration period equal to half the
    // current one, not less.
    if (z.tail(1) + dz.tail(1) < 0.5 * z.tail(1))
        lambda = -0.5 * z.tail(1) / dz.tail(1);

    // initial error
    double e_norm_zero   = error_norm(G, S, z, tmp);
    double e_norm_lambda = 0.0;

    for (auto ls_iter = 1; ls_iter <= opts.ls_maxiter; ls_iter++) {
        // calc error
        e_norm_lambda = error_norm_lambda(G, S, lambda, z, dx, tmp);

        // might exit earlier
        if (e_norm_lambda < e_norm_zero) {
            return std::make_tuple(lambda, e_norm_lambda);
        }

        lambda *= opts.ls_rho;
    }

    throw std::runtime_error("ls_maxiter reached")
}

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
        auto [lambda, e_norm] = linesearch(G, S, z, nrp.rhs(), tmp);

        // add correction
        z = z + lambda * rnp.rhs();

        // print output
        std::cout << iter << " " << e_norm << " " << lambda << "\n";

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