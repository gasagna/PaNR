#pragma once

#include "utils.hpp"

namespace PaNr {

template <typename L, typename X>
class LinOperatorExec {
private:
    LinOperator<L, X>& _op;           // actual operator
    X&                 _y;            // linearised state to be mapped
    bool               _is_last_rank; // whether to apply spatial shifts

public:
    // constructor
    LinOperatorExec(LinOperator<L, X>& op, X& _y, const bool is_last_rank)
        : _op(op)
        , _y(y)
        , _is_last_rank(is_last_rank) {}

    // actually execs the action of the operator, overwriting dest
    void execute(const X& dest) {
        dest   = _y;                                      // copy to destination
        _op._x = _op._z.head();                           // then copy initial nonlinear state
        _op._l(op._x, dest, 0, _op._z.tail(1));           // then propagate
        if (_is_last_rank)                                // only apply if needed
            static_for_enum(_op._s,
                [&](auto i, auto& s) { s(dest, _op._z.tail(i)); }) // apply all shifts
    }
};

template <typename L, typename S, typename X>
class LinOperator {
private:
    L&                   _l; // Linearised coupled operator
    S&                   _s; // Shift operator (only for last interval)
    DVector<X, NBORDER>& _z; // Nonlinear initial state (not owning, just a reference).
                             // This also contains the time and spatial shifts
    X _x;                    // temporary for the forward propagation

public:
    // constructor
    LinOperator(L& l, S& s, DVector<X, NBORDER>& z)
        : _l(l)           // linear operator
        : _s(s)           // tuple of spatial shift operators
        , _z(z)           // just a reference
        , _x(z.head()) {} // makes copy from z

    // lazy dot product
    LinOperatorExec<L, X> operator*(const X& y) {
        bool is_last_rank = _z.dinfo().rank() == _z.dinfo().comm_size() - 1;
        return { *this, y, is_last_rank }
    }

    // allow LinOperatorExec to access the provate data in here
    friend class LinOperatorExec;
};
}