 #include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(BH)]]
#include "boost/numeric/odeint.hpp"
namespace odeint = boost::numeric::odeint;

__HEADERS__;

namespace odeintr
{
  static const std::size_t N = __SYS_SIZE__;
  
  using state_type = std::array<double, N>;
  
  static state_type state;
  
  using stepper_type = odeint::__STEPPER_TYPE__;
  
  static auto stepper = __STEPPER_CONSTRUCT__;
  
  using vec_type = std::vector<double>;
  static std::array<vec_type, N> rec_x;
  static std::vector<double> rec_t;
  
  __GLOBALS__;
  
  static void
  sys(const state_type x, state_type &dxdt, const double t)
  {
    __SYS__;
  }

  static void
  obs(const state_type x, const double t)
  {
    for (int i = 0; i != N; ++i)
      rec_x[i].push_back(x[i]);
    rec_t.push_back(t);
  }
}; // namespace odeintr

// [[Rcpp::export]]
Rcpp::List __FUNCNAME___get_output()
{
  Rcpp::List out;
  out("t") = Rcpp::wrap(odeintr::rec_t);
  for (int i = 0; i != odeintr::N; ++i)
  {
    auto cnam = std::string("x") + std::to_string(i + 1);
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
  __FUNCNAME___reset_observer();
  odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state,
                             start, start + duration, step_size,
                             odeintr::obs);
  return __FUNCNAME___get_output();
}

// [[Rcpp::export]]
Rcpp::List __FUNCNAME___at(Rcpp::NumericVector init,
                           std::vector<double> times,
                           double step_size = 1.0)
{
  __FUNCNAME___set_state(init);
  __FUNCNAME___reset_observer();
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
  __FUNCNAME___reset_observer();
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
Rcpp::List
__FUNCNAME___continue(double duration, double step_size = 1.0)
{
  double start = odeintr::rec_t.back();
  odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state,
                             start, start + duration, step_size, odeintr::obs);
  return __FUNCNAME___get_output();
}

__FOOTERS__;


