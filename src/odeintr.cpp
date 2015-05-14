// Copyright Timothy H. Keitt 2015

#include "../inst/include/odeintr.h"

static std::vector<List> rec_x;
static std::vector<double> rec_t;

// using stepper_type = odeint::runge_kutta_dopri5<state_type>;
typedef odeint::runge_kutta_dopri5<state_type> stepper_type;

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

static void
init_obs(int sz)
{
  rec_x.resize(0);
  rec_t.resize(0);
  rec_x.reserve(sz);
  rec_t.reserve(sz);  
}

// [[Rcpp::export]]
List integrate_sys_const(Function derivs, Function recfun, state_type init,
                         double duration, double step_size = 1.0, double start = 0.0,
                         double atol = 1e-6, double rtol = 1e-6)
{
  Sys s(derivs); Obs o(recfun);
  init_obs(duration / step_size);
  auto stepper = odeint::make_dense_output(atol, rtol, stepper_type());
  odeint::integrate_const(stepper, s, init, start, start + duration, step_size, o);
  List out; out("t") = rec_t; out("x") = rec_x;
  return out;
}

// [[Rcpp::export]]
List integrate_sys_adapt(Function derivs, Function recfun, state_type init,
                         double duration, double step_size = 1.0, double start = 0.0,
                         double atol = 1e-6, double rtol = 1e-6)
{
  Sys s(derivs); Obs o(recfun);
  init_obs(duration / step_size);
  auto stepper = odeint::make_dense_output(atol, rtol, stepper_type());
  odeint::integrate_adaptive(stepper, s, init, start, start + duration, step_size, o);
  List out; out("t") = rec_t; out("x") = rec_x;
  return out;
}

