#include "boost/numeric/odeint.hpp"

#include <Rcpp.h>
using namespace Rcpp;

namespace odeint = boost::numeric::odeint;

struct Sys
{
  Sys(Function f) : derivs(f) {}
  void operator()(const NumericVector &x, NumericVector &dxdt, double t) const
  {
    NumericVector res = derivs(x, t);
    std::copy(res.begin(), res.end(), dxdt.begin());
  }
  Function derivs;
};

struct Obs
{
  Obs(Function f, Function g) : rec(f), comb(g) {}
  void operator()(const NumericVector &x, double t)
  {
    log = comb(log, rec(x, t));
  }
  SEXP get_log() { return wrap(log); }
  Function rec, comb;
  SEXP log;
};

//' @export
// [[Rcpp::export]]
SEXP integrate(Function derivs, Function obs, Function comb,
               NumericVector init, double from, double to, double by)
{
  Sys s(derivs);
  Obs o(obs, comb);
  odeint::integrate(s, init, from, to, by, o);
  return o.get_log();
}
