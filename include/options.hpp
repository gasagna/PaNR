#pragma once

namespace PaNr {
struct Options {
    int    maxiter        = 1.0;
    double init_tr_radius = 1.0;
    double max_tr_radius  = 1.0;
    double eta            = 0.25;
    int    gmres_maxiter  = 50;
    double gmres_res_tol  = 1e-4;
    bool   gmres_verbose  = true;
};
}