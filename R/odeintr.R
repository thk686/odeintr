#' Odeintr: Fast and Flexible Integration of Ordinary Differential Equations
#'
#' This package is a light-weight wrapper around the Boost ODEINT
#' library. It allows one to specify an ODE system in a few lines
#' of C++. This code is inserted into a template that is compiled.
#' The resulting Rcpp function will integrate the system.
#' 
#' You can also specify the model as an R function. Unlike
#' existing packages, you can also supply an observer function
#' that can return arbitrary data structures.
#' 
#' The main function is \code{compile_sys}, which takes a
#' snippet of C++ code calculating dervatives and compiles
#' an integrator function.
#' 
#' The function \code{integrate_sys} accepts an R function
#' defining the system and an observer function to record
#' the output in a data frame or list.
#' 
#' @author Timothy H. Keitt \cr \url{http://www.keittlab.org/} \cr \cr
#' 
#' Timothy H. Keitt \email{tkeitt@@gmail.com} \cr
#' 
#' @references \url{http://github.com/thk686/odeintr}, \url{http://headmyshoulder.github.io/odeint-v2/}
#' 
#' @keywords package
#' 
#' @import Rcpp
#' 
#' @useDynLib odeintr
#' 
#' @docType package
#' 
#' @name odeintr
#' @rdname odeintr_package
NULL

#' Integrate an ODE system using ODEINT
#' 
#' Numerically integrates an ODE system defined in R
#' 
#' @param dxdt a function with signature function(x, t)
#' @param init the initial conditions
#' @param duration time-span of the integration
#' @param step_size the initial step size (adjusted internally)
#' @param start the starting time
#' @param observer a function with signature function(x, t) returning values to store in output
#' 
#' @details The system will be integrated from \code{start} to \code{start + duration}. The method
#' is an error controlled 5th-order Dormand-Prince. The time step will be adjusted to within error
#' tolerances (1e-6 absolute and relative).
#' 
#' The observer can return arbitrary data in any form that can be coerced to a list. This could
#' be a single scalar value (no need to wrap the return with \code{list}!) or a list containing
#' heterogeneous types. These will be inserted into the columns of the returned data frame. If
#' the observer function returns a zero-length object (\code{NULL} or \code{list()}), then nothing
#' will be recorded. You can use the \code{t} argument to selectively sample the output.
#' 
#' @return A data frame, \code{NULL} if no samples were recorded and a very complicated
#' list-of-lists if the observer returned objects of different length.
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{compile_sys}}
#' 
#' @examples
#' \dontrun{
#' # Lotka-Volterra predator-prey equations
#' LV.sys = function(x, t)
#' {
#'    c(x[1] - 0.1 * x[1] * x[2],
#'      0.05 * x[1] * x[2] - 0.5 * x[2])
#' }
#' null_rec = function(x, t) NULL
#' system.time(integrate_sys(LV.sys, rep(1, 2), 1e3, observer = null_rec))
#' named_rec = function(x, t) c(Prey = x[1], Predator = x[2])
#' x = integrate_sys(LV.sys, rep(1, 2), 100, observer = named_rec)
#' plot(x[, 2:3], type = "l", lwd = 3, col = "steelblue")
#' Sys.sleep(0.5)
#' 
#' # Lorenz model from odeint examples
#' Lorenz.sys = function(x, t)
#' {
#'  c(10 * (x[2] - x[1]),
#'    28 * x[1] - x[2] - x[1] * x[3],
#'    -8/3 * x[3] + x[1] * x[2])
#' }
#' system.time(integrate_sys(Lorenz.sys, rep(1, 3), 1e2, obs = null_rec))
#' x = integrate_sys(Lorenz.sys, rep(1, 3), 100)
#' plot(x[, c(2, 4)], type = 'l', col = "steelblue")
#' }
#' @export
integrate_sys = function(sys, init, duration,
                         step_size = 1, start = 0,
                         observer = function(x, t) x)
{
  res = integrate_sys_(sys, observer, init, duration, step_size, start)
  if (length(res[[1]]) == 0) return(NULL)
  x = res[[2]]; out = list()
  if (any(diff(sapply(x, length)) != 0)) return(res)
  n = length(x[[1]]); length(out) = n + 1
  for (i in 1:length(x)) for (j in 1:n)
      out[[j + 1]][i] = x[[i]][[j]]
  out[[1]] = res[[1]]
  xnames = names(x[[1]])
  if (is.null(xnames) || length(xnames) != length(x[[1]]))
    xnames = paste0("x", 1:length(x[[1]]))
  names(out) = c("t", xnames)
  as.data.frame(out)
}
  
#' Compile ODE system
#' 
#' Generates an integrator using Rcpp
#' 
#' @param name the name of the generated integration function
#' @param sys a string containing C++ expressions
#' @param globals a string with global C++ declarations
#' @param headers code to appear before the \code{odeintr} namespace
#' @param footers code to appear after the \code{odeintr} namespace
#' @param sys_dim length of the state vector
#' @param compile if false, just return the code
#' @param ... passed to \code{\link{sourceCpp}}
#' 
#' @details C++ code is generated and compiled with
#' \code{\link{sourceCpp}}. The returned function will
#' integrate the system starting from a provided initial
#' condition and initial time to a specified final time.
#' An attempt is made to get the length of the state vector
#' from the system definition. If this fails, the code will
#' likely crash your R session. It is safer to specify
#' \code{sys_dim} directly.
#' 
#' The signature of \code{name} is:
#' 
#' \code{name(initial_state, duration, step_size, start)}
#' 
#' where \code{duration} is the integration time, \code{step_size}
#' is the time step (default = 1.0) and \code{start} is the
#' starting time (default = 0.0). The
#' return value from this function is a data frame with time
#' in the first column and state in the remaining columns.
#' 
#' A second function, \code{name_no_record}, is also generated.
#' This function does not take an observer and will return the
#' approximate solution at the end of the integration. Use this
#' for benchmarking or discarding transients (result can be used
#' as the initial condition for a subsequent call).
#' 
#' The integration method is a 5th order Dormand-Prince algorithm
#' with error tolerance (absolute and relative) of 1e-6. The
#' step size is adaptive. The given step size is tried and a
#' smaller or larger step size will be used depending on the
#' esimtated integration error. The time and system state are
#' recorded after every step. The time increments will generally
#' be variable.
#' 
#' The C++ expressions must index a state array of length
#' \code{sys_dim}. The state array is \code{x} and the
#' derivatives are \code{dxdt}. The first state value is
#' \code{x[0]} and the first derivative is \code{dxdt[0]}.
#' 
#' In the case you use bare \code{dxdt} and \code{x}, an
#' attempt will be made to append \code{[0]} to these
#' variables. This can fail, so do not rely on it. It
#' requires that \code{dxdt} and \code{x} are adjacent
#' to a "non-word" character as defined by perl regular
#' expressions. This will also fail if you set \code{sys_dim}
#' to a positive value.
#' 
#' The \code{globals} string can be arbitrary C++ code. You
#' can set global named parameter values here. Note that
#' these will be defined within the \code{odeintr} namespace.
#' 
#' You can insert arbitrary code compiled outside the \code{odeintr}
#' names space using \code{headers} and \code{footers}. This code
#' can be anything compatible with Rcpp. You could for example
#' define exported Rcpp functions that set simulation paramters.
#' \code{headers} is inserted right after the Rcpp and ODEINT
#' includes. \code{footers} is inserted at the end of the
#' code.
#' 
#' @note The c++11 plugin is enabled.
#' 
#' Also, despite using an array as the state,
#' even when it is length one, there is no performance penalty as the length
#' of the state vector is known at compile time. Your compiler will optimize
#' out any pointer dereferencing. You may get a large increase in performance
#' by enabling function inlining and loop unrolling. For GCC, that is "-O3". Chances
#' are that your R distribution uses "-O2". You can override this by copying each
#' line containing "-O2" from your R distribution's Makeconf to .R/Makevars in
#' your home directory replacing each instance of "-O2" with "-O3". You might
#' be able to do this with environment variables. I have not been able to
#' verify whether that works. You might also try \code{\link{set_optimization}}.
#'  
#' @return the C++ code invisibly
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{set_optimization}}, \code{\link{integrate_sys}}
#' 
#' @examples
#' \dontrun{
#' # Logistic growth
#' compile_sys("logistic", "dxdt = x * (1 - x)")
#' plot(logistic(0.001, 15), type = "l", lwd = 2, col = "steelblue")
#' Sys.sleep(0.5)
#' 
#' # Lotka-Volterra predator-prey equations
#' LV.sys = '
#'   dxdt[0] = x[0] - 0.1 * x[0] * x[1];
#'   dxdt[1] = 0.05 * x[0] * x[1] - 0.5 * x[1];
#' ' # LV.sys
#' compile_sys("preypred", LV.sys)
#' system.time(preypred_no_record(rep(1, 2), 1e6))
#' x = preypred(rep(1, 2), 100)
#' plot(x[, 2:3], type = "l", xlab = "Prey", ylab = "Predator", col = "steelblue")
#' Sys.sleep(0.5)
#' 
#' # Lorenz model from odeint examples
#' Lorenz.globals = '
#'   const double sigma_ = 10.0;
#'   const double R_ = 28.0;
#'   const double b_ = 8.0 / 3.0;
#' ' # Lorenz.globals
#' Lorenz.sys = '
#'   dxdt[0] = sigma_ * (x[1] - x[0]);
#'   dxdt[1] = R_ * x[0] - x[1] - x[0] * x[2];
#'   dxdt[2] = -b_ * x[2] + x[0] * x[1];
#' ' # Lorenz.sys
#' compile_sys("lorenz", Lorenz.sys, Lorenz.globals)
#' system.time(lorenz_no_record(rep(1, 3), 1e5))
#' x = lorenz(rep(1, 3), 100)
#' plot(x[, c(2, 4)], type = 'l', col = "steelblue")
#' }
#' @export
compile_sys = function(name, sys,
                       method = "rk5_i",
                       sys_dim = -1L,
                       atol = 1e-6,
                       rtol = 1e-6,
                       globals = "", 
                       headers = "",
                       footers = "",
                       compile = TRUE,
                       ...)
{
  stepper_constr = make_stepper_constr(method, atol, rtol)
  stepper = make_stepper_type(method)
  if (ceiling(sys_dim) < 1) sys_dim = get_sys_dim(sys)
  if (is.na(sys_dim))
  {
    sys = gsub("\\Wdxdt|dxdt\\W", "dxdt[0]", sys, perl = TRUE)
    sys = gsub("\\Wx|x\\W", "x[0]", sys, perl = TRUE)
    sys_dim = 1L
  }
  fn = system.file(file.path("templates", "rcpp_template.cpp"),
                   package = "odeintr", mustWork = TRUE)
  con = file(fn); code = readLines(con); close(con)
  code = gsub("__FUNCNAME__", name, code)
  code = gsub("__STEPPER_TYPE__", stepper, code)
  code = gsub("__STEPPER_CONSTRUCT__", stepper_constr, code)
  code = sub("__SYS_SIZE__", ceiling(sys_dim), code)
  code = sub("__GLOBALS__", globals, code)
  code = sub("__SYS__", sys, code)
  code = sub("__HEADERS__", headers, code)
  code = sub("__FOOTERS__", footers, code)
  code = paste0(code, collapse = "\n")
  if (compile)
    tryCatch(Rcpp::sourceCpp(code = code, ...),
             error = function(e) e)
  return(invisible(code))
}

make_stepper_constr = function(method, atol, rtol)
{
  stepper_constr = "stepper_type()"
  if (grepl("_a$", method))
    stepper_constr = paste0("odeint::make_controlled(", atol, ", ", rtol, ", stepper_type())")
  if (grepl("_i$", method))
    stepper_constr = paste0("odeint::make_dense_output(", atol, ", ", rtol, ", stepper_type())") 
  return(stepper_constr)
}

make_stepper_type = function(stepper)
{
  if (grepl("ab|am|abm", stepper))
  {
    steps = as.integer(sub("ab|am|abm([0-9]+)", "\\1", stepper))
    stepper = sub("(ab|am|abm)[0-9]+", "\\1", stepper)
  }
  switch(stepper,
         euler = "euler<state_type>",
         rk4 = "runge_kutta4<state_type>",
         rk54 = "runge_kutta_cash_karp54<state_type>",
         rk5 = "runge_kutta_dopri5<state_type>",
         rk78 = "runge_kutta_fehlberg78<state_type>",
         ab = paste0("adams_bashforth<", steps, ", state_type>"),
         am = paste0("adams_moulton<", steps, ", state_type>"),
         abm = paste0("adams_bashforth_moulton<", steps, ", state_type>"),
         bs = "bulirsch_stoer<state_type>",
         bsd = "bulirsch_stoer_dense_out<state_type>",
         paste0(stepper, "<state_type>"))
}

get_sys_dim = function(x)
{
  matches = gregexpr("dxdt\\[\\d\\]", x, perl = TRUE)
  lens = attr(matches[[1]], "match.length") - 7L
  starts = unlist(matches) + 5L
  indices = rep(NA, length(starts))
  for (i in seq(along = indices))
    indices[i] = substr(x, starts[i], starts[i] + lens)
  return(max(as.integer(indices)) + 1L)
}

array_sys_template = function()
{
  '
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
  Rcpp::List get_output()
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
  Rcpp::List __FUNCNAME__(Rcpp::NumericVector init,
                          double duration,
                          double step_size = 1.0,
                          double start = 0.0)
  {
    if (init.size() != odeintr::N)
      Rcpp::stop("Invalid initial state");
    odeintr::state_type inival;
    std::copy(init.begin(), init.end(), inival.begin());
    for (auto &i : odeintr::rec_x) i.resize(0);
    odeintr::rec_t.resize(0);
    odeint::integrate(odeintr::sys, inival,
                      start, start + duration, step_size,
                      odeintr::obs);
    return get_output();
  }

  // [[Rcpp::export]]
  std::vector<double>
  __FUNCNAME___no_record(Rcpp::NumericVector init,
                         double duration,
                         double step_size = 1.0,
                         double start = 0.0)
  {
    if (init.size() != odeintr::N) Rcpp::stop("Invalid initial state");
    odeintr::state_type inival;
    std::copy(init.begin(), init.end(), inival.begin());
    odeint::integrate(odeintr::sys, inival, start, start + duration, step_size);
    return std::vector<double>(inival.begin(), inival.end());
  }
  
  __FOOTERS__;
  '
}

disable_asserts = function(makevars)
{
  con = pipe(paste("R CMD config CPPFLAGS"))
  flags = readLines(con); close(con)
  if (!is.finite(pmatch("-DNDEBUG", flags)))
    flags = paste(flags, "-DNDEBUG")
  flags = gsub("^\\s+|\\s+$", "", flags, perl = TRUE)
  cat(paste0("CPPFLAGS=", flags, "\n"),
      file = makevars, append = TRUE)
}

substitute_opt_level = function(flags, level, omit.debug)
{
  flags = gsub("-O\\d+", paste0("-O", level), flags, perl = TRUE)
  if (omit.debug) flags = gsub("\\s*-g\\s*", "", flags, perl = TRUE)
  flags = gsub("^\\s+|\\s+$", "", flags, perl = TRUE)
  return(flags)
}

process_flags = function(name, level, omit.debug)
{
  con = pipe(paste("R CMD config", name))
  flags = readLines(con); close(con)
  flags = substitute_opt_level(flags, level, omit.debug)
  paste0(name, "=", flags)
}

.opt.env.vars = c("CFLAGS", "FFLAGS", "FCFLAGS", "CXXFLAGS", "CXX1XFLAGS")

#' Set compiler optimization
#' 
#' Write a user Makevars with updated optimization level
#' 
#' @param level the compiler optimization level (-O<level>)
#' @param omit.debug if true, remove "-g" from flags
#' @param disable.asserts if true, set NDEBUG define
#' @param overwrite if true, overwrite existing Makevars file
#' 
#' @details This function will change the optimization flags used
#' when compiling code. It will write the file "Makevars" to the
#' ".R" directory in your "$HOME" directory. These setting will
#' effect all subsequent compiles, even package installation,
#' unless you remove or edit the "Makevars" file.
#' 
#' This function assumes that your compiler uses "-O" to
#' indicate optimization level and "-g" to indicate that
#' the compiler should issue debugging symbols.
#' 
#' Please let me know if this function fails for your
#' compiler. (Submit and issue on GitHub.)
#' 
#' @note Don't go overboard. Levels greater than 3 can be
#' hazardous to numerical accuracy. Some packages will not
#' compile or will give inaccurate results for levels above 2.
#' 
#' A very similar function exists in the RStan package.
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{compile_sys}}
#' 
#' @rdname set-opt
#' @export
set_optimization = function(level = 3,
                            omit.debug = FALSE,
                            disable.asserts = FALSE,
                            overwrite = FALSE)
{
  user_dir = file.path(Sys.getenv("HOME"), ".R")
  if (!file.exists(user_dir)) dir.create(user_dir)
  makevars = file.path(user_dir, "Makevars")
  if (overwrite && file.exists(makevars)) unlink(makevars)
  if (file.exists(makevars))
    stop("User Makevars file exists; use overwrite = TRUE")
  lapply(.opt.env.vars, function(x)
    {
      cat(process_flags(x, level, omit.debug), "\n",
          file = makevars, append = TRUE)
    })
  if (disable.asserts) disable_asserts(makevars)
  invisible(.opt.env.vars)
}

#' @rdname set-opt
#' @export
show.Makevars = function()
{
  user_dir = file.path(Sys.getenv("HOME"), ".R")
  makevars = file.path(user_dir, "Makevars")
  mv = readLines(makevars)
  if (file.exists(makevars))
    for (l in mv) cat(l, "\n")
  return(invisible(mv))
}

#' @rdname set-opt
#' @export
rm.Makevars = function()
{
  user_dir = file.path(Sys.getenv("HOME"), ".R")
  makevars = file.path(user_dir, "Makevars")
  if (file.exists(makevars))
  {
    mv = readLines(makevars)
    unlink(makevars)
  }
  return(invisible(mv))
}

