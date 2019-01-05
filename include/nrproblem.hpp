#pragma once

#include "linop.hpp"
#include "utils.hpp"

namespace PaNr {

template <typename GT, typename LT, typename DGT, typename DST>
class NRProblem {
private:
    GT&                                     _g;
    LT&                                     _l;
    ST&                                     _s;
    dG&                                     _dg;
    dS&                                     _ds;
    LinOperator<LT, X>                      _oper;
    DMatrix<LinOperator<LT, X>, X, NBORDER> _lhs;
    DVector<X, NBORDER>                     _rhs;

public:
    NRProblem(GT& g, LT& l, ST& s, dG& dg, dS& ds, DVector<X, NBORDER>& z0)
        : _g(g)
        , _l(l)
        , _s(s)
        , _dg(g)            // this must be a tuple of length NBORDER
        , _oper(l, s, z0)   // _oper will hold a reference to z0
        , _tmp(z0.head())   // These are three temporaries used for the
        , _xT(z0.head())    // matrix vector product (see update)
        , _dxTdT(z0.head()) //
        , _rhs(z0)          // a copy for the right hand side
        , _lhs(z0.dinfo().comm(), _oper, DMatrixBandType::UPPER, z0) {

        static_assert(_dg.size() == NBORDER);
        static_assert(_s.size() == NBORDER - 1);
    }

    auto& rhs() { return _rhs; }
    auto& lhs() { return _lhs; }

    void update() {
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        // We first need to understand where we are. So get number of processors
        // across which z0 is distributed and the local processor rank
        int N = _z0.dinfo().size();
        int i = _z0.dinfo().rank();

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        // Then, make sure we get the right values for the tail. This is
        // used in a few place, for instance by _oper, to get the time for
        // the forward map or below to propagate states forward and
        // calculate various quantities at the end of the segment. The
        // operators on the diagonal are fed directly by _z0, so we do not need
        // to update them.
        _z0.bc_tail();

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        // Fill the lower bordering vectors with the time and shift derivatives
        // of the initial state in the first slot, while the rest is set to zero
        if (i == 0) {
            // we need this because we've got a tuple of functions
            static_for_enum(_dg,
                            [&](auto  j,
                                auto& dg) { dg(_rhs.dborders(j).head(),
                                               _z0.head()); });
        } else {
            for (auto j = 0; j != NBORDER; j++)
                nrp.rhs().dborders(j).head() = 0;
        }

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        // now update final state by first making a copy of the local
        // initial state and then propagating it forward by T/N
        _xT = _z0.head();
        _g(_xT, 0, _z0.tail(1) / N);

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        // obtain finite difference approximation of
        // derivative of the flow operator g with a second integration
        _tmp = _z0.head();
        _g(_tmp, 0, _z0.tail(1) / N + _epsilon_T);
        _rhs.rborders(0).head() = (_tmp - _xT) / _epsilon_T;

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        // apply shifts to the last segment
        if (i == N - 1) {
            static_for_enum(_s,
                            [&](auto  j,
                                auto& s) { _s(_xT,    _z0.tail(j+1)); 
                                           _s(_dxTdT, _z0.tail(j+1)); });
        }
    }
};
}