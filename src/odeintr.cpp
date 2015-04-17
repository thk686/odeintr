#include "../inst/include/odeintr.h"

// Globals required because observer is passed by value
static std::vector<List> rec_x;
static std::vector<double> rec_t;

struct Sys
{
  Sys(Function f) : derivs(f) {}
  void operator()(const state_type &x, state_type &dxdt, double t) const
  {
    auto d = as<state_type>(derivs(x, t));
    if (d.size() != dxdt.size()) stop("Invalid dimensions");
    std::copy(d.begin(), d.end(), dxdt.begin());
  }
  Function derivs;
};

struct Obs
{
  Obs(Function f) : recf(f) {}
  void operator()(const state_type x, const double t)
  {
    List rec = recf(x, t);
    if (rec.length() != 0)
    {
      rec_x.push_back(rec);
      rec_t.push_back(t);
    }
  }
  Function recf;
};

// [[Rcpp::export]]
List integrate_sys_(Function derivs, Function recfun, state_type init,
                    double duration, double step_size = 1.0,
                    double start = 0.0)
{
  rec_x.resize(0);
  rec_t.resize(0);
  Sys s(derivs);
  Obs o(recfun);
  odeint::integrate(s, init, start, start + duration, step_size, o);
  List out;
  out("t") = rec_t;
  out("x") = rec_x;
  return out;
}

