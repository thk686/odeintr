make_sys = function(sys)
{
  res = "
    struct Sys
    {
      void operator()(const std::vector<double> &x, std::vector<double> &dxdt, double t) const
      {
        __the_system__;
      }
    };
  "
  res = sub("__the_system__", sys, res)
  return(res)
}

make_obs = function(state_type, rec, comb, get = "wrap(state)")
{
  res = "
    struct Obs
    {
      void operator()(const std::vector<double> &x, double t)
      {
        state = comb(state, rec(x, t));
      }
      __the_recorder__;
      __the_combinator__;
      SEXP get() { __the_getter__; }
      __the_state_type__ state;
    };
  "
  res = gsub("__the_state_type__", state_type, res)
  res = sub("__the_combinator__", comb, res)
  res = sub("__the_getter__", get, res)
  return(res)
}

matrix_observer = function(state_size)
{
  res = "
    struct Obs
    {
      Obs() : state(__state_size__ + 1) {}
      void operator()(const std::vector<double> &x, double t)
      {
        if (x.size() != __state_size__) stop(\"Wrong state size\");
        for (int i = 0; i != __state_size__; ++i)
          state[i + 1].push_back(x[i]);
        state[0].push_back(t);
      }
      SEXP get() { Rcpp::wrap(state); }
      std::vector<std::vector<double> > state;
    };
  "
  res = gsub("__state_size__", state_size, res)
  return(res)  
}

make_code = function(sys, obs, fname = "ode")
{
  res = "
    #include <Rcpp.h>
    using namespace Rcpp;

    // [[Rcpp::depends(BH)]]
    #include \"boost/numeric/odeint.hpp\"
    namespace odeint = boost::numeric::odeint;

    __the_system__

    __the_observer__

    //' @export
    // [[Rcpp::export]]
    SEXP __the_function_name__(std::vector<double> init, double t0, double tn, double ts)
    {
      Sys s;
      Obs o;
      odeint::integrate(s, init, t0, tn, ts, o);
      return Rcpp::wrap(o.get());
    }
  "
  res = sub("__the_system__", sys, res)
  res = sub("__the_observer__", obs, res)
  res = sub("__the_function_name__", fname, res)
  return(res)
}

compile_sys = function(sys_dim = 1,
                       sys = "dxdt[0] = x[0] * (1 - x[0])",
                       funcname = "run_it",
                       ...)
{
  code = array_sys_template()
  code = sub("__FUNCNAME__", funcname, code)
  code = sub("__SYS_SIZE__", sys_dim, code)
  code = sub("__SYS__", sys, code)
  Rcpp::sourceCpp(code = code, ...)
  invisible(code)
}

array_sys_template = function()
{
  '
    #include <Rcpp.h>
    using namespace Rcpp;
    // [[Rcpp::plugins(cpp11)]]
    
    // [[Rcpp::depends(BH)]]
    #include "boost/numeric/odeint.hpp"
    namespace odeint = boost::numeric::odeint;
    
    static const std::size_t N = __SYS_SIZE__;
    
    using state_type = std::array<double, N>;
    using vec_type = std::vector<double>;
    
    static std::array<vec_type, N> m_x;
    static std::vector<double> m_t;
    
    static void
    sys(const state_type x, state_type &dxdt, const double t)
    {
      __SYS__;
    }
      
    static void
    obs(const state_type x, const double t)
    {
      for (int i = 0; i != N; ++i)
        m_x[i].push_back(x[i]);
      m_t.push_back(t);
    }
      
    // [[Rcpp::export]]
    List __FUNCNAME__(NumericVector init,
                      double from, double to,
                      double by)
    {
      if (init.size() != N)
        stop("Invalid initial state");
      state_type inival;
      for (int i = 0; i != N; ++i) inival[i] = init[i];
      odeint::integrate(sys, inival, from, to, by, obs);
      List out;
      out("t") = wrap(m_t);
      for (int i = 0; i != N; ++i)
      {
        auto cnam = std::string("x") + std::to_string(i + 1);
        out(cnam) = wrap(m_x[i]);
      }
      out.attr("row.names") = IntegerVector::create(NA_INTEGER, -m_t.size());
      out.attr("class") = "data.frame";
      return out;
    }
  '
}
