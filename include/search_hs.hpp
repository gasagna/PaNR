
/* main search function */
template <typename GT,
          typename LT,
          typename FT,
          typename ST,
          typename DST,
          typename X,
          std::size_t NBORDER>
void _search_hookstep(GT&&                 G, // nonlinear propagator
                      LT&&                 L, // linearised propagator
                      ST&&                 S, // tuple of shift operators (can be empty tuple)
                      DT&&                 D, //
                      DVector<X, NBORDER>& z0,
                      Options              opt) {

    // define object for solving the NR iterations, coupling it
    // to z0. Every change in z0 is seen by nrp.
    auto nrp = NRProblem(G, L, S, D, z0);

    // init trust region radius and other parameters
    double tr_radius = opts.init_tr_radius;
    int    iter      = 0;
    double rho       = 1.0;

    // main looping
    while (true) {
        // update the linear problem with the current solution
        nrp.update();

        // solve linear problem in place, overwriting the rhs
        solve_gmres(nrp.lhs(),
                    nrp.rhs(),
                    tr_radius, // hookstep!
                    opts.gmres_rtol,
                    opts.gmres_maxiter,
                    opts.gmres_verbose);

        // calc ratio of actual and predicted decrease
        rho = actual_over_predicted();

        // update trust region radius
        if (rho < 0.25)
            tr_radius *= 0.25;
        else if (rho > 0.75 && hits_boundary)
            tr_radius = std::min(2 * tr_radius, opts.max_tr_radius);

        // update current solution
        if (rho > opts.eta)
            z0 += rnp.rhs();

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