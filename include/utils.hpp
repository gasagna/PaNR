#pragma once
#include <tuple>
#include <utility>

namespace PaNr {

//////////////////////////////////////////////////////////////////////////
// Utilities for static iteration over tuples. See
// https://codereview.stackexchange.com/questions/173564/
// This uses c++17 fold expressions

//////////////////////////////////////////////////////////////////////////
//
// Allow enumeration over a tuple, like in Python
//
// for i, el in enumerate(tup):
//      f(i, el)
//
// helper
template <typename TUP, typename FUN, std::size_t... I>
constexpr void _static_for_enum(TUP& tup, FUN&& f, std::index_sequence<I...>) {
    (f(std::integral_constant<std::size_t, I>(), std::get<I>(tup)), ...);
}

// main function
template <typename... T, typename FUN>
constexpr void static_for_enum(std::tuple<T...>& tup, FUN&& f) {
    _static_for_enum(tup,
                     std::forward<FUN>(f),
                     std::make_index_sequence<sizeof...(T)>());
}

// simple overload for non tuples
template <typename T, typename FUN>
constexpr void static_for_enum(T& el, FUN&& f) { f(0, el); }

//////////////////////////////////////////////////////////////////////////
//
// Allow iteration over a tuple, like in Python
// for el in tup:
//      f(el)
//
// helper
template <typename TUP, typename FUN, std::size_t... I>
constexpr void _static_for(TUP& tup, FUN&& f, std::index_sequence<I...>) {
    (f(std::get<I>(tup)), ...);
}

// main function
template <typename... T, typename FUN>
constexpr void static_for(std::tuple<T...>& tup, FUN&& f) {
    _static_for(tup,
                std::forward<FUN>(f),
                std::make_index_sequence<sizeof...(T)>());
}

// simple overload for non tuples
template <typename T, typename FUN>
constexpr void static_for(T& el, FUN&& f) { f(el); }
}