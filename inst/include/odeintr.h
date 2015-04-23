// Copyright Timothy H. Keitt 2015

#ifndef __ODEINTR_H__
#define __ODEINTR_H__

#include "boost/numeric/odeint.hpp"

#include <Rcpp.h>
using namespace Rcpp;

namespace odeint = boost::numeric::odeint;

using state_type = std::vector<double>;

#endif // __ODEINTR_H__
