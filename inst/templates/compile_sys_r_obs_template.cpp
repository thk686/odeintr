// Copyright Timothy H. Keitt 2015
// See license for odeintr package

// [[Rcpp::depends(odeintr)]]

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(BH)]]
#include "boost/numeric/odeint.hpp"
namespace odeint = boost::numeric::odeint;

__HEADERS__;

namespace odeintr
{
  static const std::size_t N = __SYS_SIZE__;
  
  typedef std::vector<double> state_type;
  
  static state_type state(N);
  
  typedef odeint::__STEPPER_TYPE__ stepper_type;
  
  static auto stepper = __STEPPER_CONSTRUCT__;
  
  typedef Rcpp::List vec_type;
  typedef std::vector<vec_type> rec_type;
  static rec_type rec_x;
  static std::vector<double> rec_t;
  
  static Rcpp::Function recf("c"), proc("c");
  
  __GLOBALS__;
  
  #include "utils.h"
  
  static void
  sys(const state_type x, state_type &dxdt, const double t)
  {
    __SYS__;
  }

  static void
  obs(const state_type x, const double t)
  {
    Rcpp::List rec = recf(x, t);
    if (rec.length() != 0)
    {
      rec_x.push_back(rec);
      rec_t.push_back(t);
    }
  }
  
}; // namespace odeintr

static void
reserve(int n)
{
  odeintr::rec_t.reserve(n);
  odeintr::rec_x.reserve(n);
}

// [[Rcpp::export]]
Rcpp::List __FUNCNAME___get_output()
{
  Rcpp::List out;
  out("Time") = Rcpp::wrap(odeintr::rec_t);
  out("X") = Rcpp::wrap(odeintr::rec_x);
  return odeintr::proc(out);
};

// [[Rcpp::export]]
void __FUNCNAME___set_state(Rcpp::NumericVector new_state)
{
  if (new_state.size() != odeintr::N)
    Rcpp::stop("Invalid initial state");
  std::copy(new_state.begin(),
            new_state.end(),
            odeintr::state.begin());
}

// [[Rcpp::export]]
std::vector<double>
__FUNCNAME___get_state()
{
  return odeintr::state;
}

// [[Rcpp::export]]
void __FUNCNAME___reset_observer()
{
  odeintr::rec_x.resize(0);
  odeintr::rec_t.resize(0);  
}

// [[Rcpp::export]]
Rcpp::List __FUNCNAME___adap(Rcpp::NumericVector init,
                             double duration,
                             double step_size = 1.0,
                             double start = 0.0)
{
  __FUNCNAME___set_state(init);
  __FUNCNAME___reset_observer(); reserve(duration / step_size);
  odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state,
                             start, start + duration, step_size,
                             odeintr::obs);
  return __FUNCNAME___get_output();
}

// [[Rcpp::export]]
Rcpp::List __FUNCNAME___at(Rcpp::NumericVector init,
                           std::vector<double> times,
                           double step_size = 1.0,
                           double start = 0.0)
{
  __FUNCNAME___set_state(init);
  __FUNCNAME___reset_observer(); reserve(times.size());
  odeint::integrate_const(odeintr::stepper, odeintr::sys, odeintr::state,
                          start, times[0], step_size);
  odeint::integrate_times(odeintr::stepper, odeintr::sys, odeintr::state,
                          times.begin(), times.end(), step_size, odeintr::obs);
  return __FUNCNAME___get_output();
}

// [[Rcpp::export]]
Rcpp::List
__FUNCNAME___continue_at(std::vector<double> times, double step_size = 1.0)
{
  double start = odeintr::rec_t.back();
  __FUNCNAME___reset_observer(); reserve(odeintr::rec_t.size() + times.size());
  odeint::integrate_const(odeintr::stepper, odeintr::sys, odeintr::state,
                          start, times[0], step_size);
  odeint::integrate_times(odeintr::stepper, odeintr::sys, odeintr::state,
                          times.begin(), times.end(), step_size, odeintr::obs);
  return __FUNCNAME___get_output();
}

// [[Rcpp::export]]
Rcpp::List __FUNCNAME__(Rcpp::NumericVector init,
                       double duration,
                       double step_size = 1.0,
                       double start = 0.0)
{
  __FUNCNAME___set_state(init);
  __FUNCNAME___reset_observer(); reserve(duration / step_size);
  odeint::integrate_const(odeintr::stepper, odeintr::sys, odeintr::state,
                          start, start + duration, step_size,
                          odeintr::obs);
  return __FUNCNAME___get_output();
}

// [[Rcpp::export]]
std::vector<double>
__FUNCNAME___no_record(Rcpp::NumericVector init,
                       double duration,
                       double step_size = 1.0,
                       double start = 0.0)
{
  __FUNCNAME___set_state(init);
  odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state,
                             start, start + duration, step_size);
  return __FUNCNAME___get_state();
}

// [[Rcpp::export]]
Rcpp::Function
__FUNCNAME___get_observer()
{
  return odeintr::recf;
}

// [[Rcpp::export]]
void
__FUNCNAME___set_observer(Rcpp::Function f)
{
  odeintr::recf = f;
}

// [[Rcpp::export]]
Rcpp::Function
__FUNCNAME___get_output_processor()
{
  return odeintr::proc;
}

// [[Rcpp::export]]
void
__FUNCNAME___set_output_processor(Rcpp::Function f)
{
  odeintr::proc = f;
}

__FOOTERS__;




