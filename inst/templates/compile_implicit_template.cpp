// Copyright Timothy H. Keitt 2015
// See license for odeintr package

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(BH)]]
#include "boost/numeric/odeint.hpp"
namespace odeint = boost::numeric::odeint;

#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
namespace ublas = boost::numeric::ublas;

__HEADERS__;

namespace odeintr
{
  typedef ublas::vector<double> sys_vec;
  typedef ublas::matrix<double> sys_mat;

  static const std::size_t N = __SYS_SIZE__;
  
  typedef sys_vec state_type;
  
  static state_type state(N);
  
  typedef odeint::rosenbrock4<double> stepper_type;
  
  static const double atol = __ATOL__,
                      rtol = __RTOL__;
  
  static auto stepper = odeint::make_dense_output<stepper_type>(atol, rtol);
  
  typedef std::vector<double> vec_type;
  static std::vector<vec_type> rec_x(N);
  static vec_type rec_t;
  
  __GLOBALS__;
  
  struct stiff_system
  {
      void operator()(const sys_vec &x, sys_vec &dxdt , double t)
      {
        __SYS__;
      }
  };
  
  struct stiff_system_jacobi
  {
      void operator()(const sys_vec &x, sys_mat &J,
                      const double t, sys_vec &dfdt)
      {
        __JACOBIAN__;
      }
  };

  auto sys = std::make_pair(stiff_system(), stiff_system_jacobi());

  static void
  obs(const state_type x, const double t)
  {
    for (int i = 0; i != N; ++i)
      rec_x[i].push_back(x[i]);
    rec_t.push_back(t);
  }
  
}; // namespace odeintr

static void
reserve(odeintr::vec_type::size_type n)
{
  odeintr::rec_t.reserve(n);
  for (auto &i : odeintr::rec_x) i.reserve(n);
}

// [[Rcpp::export]]
Rcpp::List __FUNCNAME___get_output()
{
  Rcpp::List out;
  out("Time") = Rcpp::wrap(odeintr::rec_t);
  for (int i = 0; i != odeintr::N; ++i)
  {
    auto cnam = std::string("X") + std::to_string(i + 1);
    out(cnam) = Rcpp::wrap(odeintr::rec_x[i]);
  }
  out.attr("class") = "data.frame";
  int rows_out = odeintr::rec_t.size();
  auto rn = Rcpp::IntegerVector::create(NA_INTEGER, -rows_out);
  out.attr("row.names") = rn;
  return out;
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
  return std::vector<double>(odeintr::state.begin(), odeintr::state.end());
}

// [[Rcpp::export]]
void __FUNCNAME___reset_observer()
{
  for (auto &i : odeintr::rec_x) i.resize(0);
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

__FOOTERS__;




