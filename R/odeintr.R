# Copyright Timothy H. Keitt 2015

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
#' @param sys a function with signature function(x, t)
#' @param init the initial conditions
#' @param duration time-span of the integration
#' @param step_size the initial step size (adjusted internally)
#' @param start the starting time
#' @param adaptive_obs if true, call observer after each adaptive step
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
#' x = integrate_sys(LV.sys, rep(1, 2), 100, 0.01, observer = named_rec)
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
#' x = integrate_sys(Lorenz.sys, rep(1, 3), 100, 0.01)
#' plot(x[, c(2, 4)], type = 'l', col = "steelblue")
#' }
#' @export
integrate_sys = function(sys, init, duration,
                         step_size = 1, start = 0,
                         adaptive_obs = FALSE,
                         observer = function(x, t) x)
{
  res = if (adaptive_obs)
    integrate_sys_adapt(sys, observer, init, duration, step_size, start)
  else
    integrate_sys_const(sys, observer, init, duration, step_size, start)
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
#' @param pars a named vector of numbers or a vector of names or number of parameters
#' @param const declare parameters const if true
#' @param method a method string (see Details)
#' @param sys_dim length of the state vector
#' @param atol absolute tolerance if using adaptive step size
#' @param rtol relative tolerance if using adaptive step size
#' @param globals a string with global C++ declarations
#' @param headers code to appear before the \code{odeintr} namespace
#' @param footers code to appear after the \code{odeintr} namespace
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
#' The C++ expressions must index a state array of length
#' \code{sys_dim}. The state array is \code{x} and the
#' derivatives are \code{dxdt}. The first state value is
#' \code{x[0]} and the first derivative is \code{dxdt[0]}.
#' 
#' In the case you use bare \code{dxdt} and \code{x}, an
#' attempt will be made to append \code{[0]} to these
#' variables. This can fail, so do not rely on it. 
#' This will also fail if you set \code{sys_dim}
#' to a positive value.
#' 
#' The \code{globals} string can be arbitrary C++ code. You
#' can set global named parameter values here. Note that
#' these will be defined within the \code{odeintr} namespace.
#' 
#' If you supply the \code{pars} argument, these parameters
#' will be compiled into the code. There are three options:
#' 1) if \code{pars} is a single number, then you can access
#' a vector of parameters named \code{pars} of the specified
#' length; 2) if \code{pars} is a character vectors, then a
#' parameter will be defined for each; and 3) if the character
#' vector is named, then the names will be used for the
#' parameter names and the associated values will be used
#' as defaults. If you specify \code{const = TRUE}, these
#' named parameters will be declared const. Otherwise
#' parameter getter/setter functions will be defined.
#' 
#' You can insert arbitrary code outside the \code{odeintr}
#' names space using \code{headers} and \code{footers}. This code
#' can be anything compatible with Rcpp. You could for example
#' define exported Rcpp functions that set simulation paramters.
#' \code{headers} is inserted right after the Rcpp and ODEINT
#' includes. \code{footers} is inserted at the end of the
#' code.
#' 
#' The following methods can be used:
#' 
#' \tabular{lll}{
#' Code \tab Stepper \tab Type \cr
#' \code{euler} \tab \code{euler} \tab Interpolating \cr
#' \code{rk4} \tab \code{runge_kutta4} \tab Regular \cr
#' \code{rk54} \tab \code{runge_kutta_cash_karp54} \tab Regular \cr
#' \code{rk54_a} \tab \code{runge_kutta_cash_karp54} \tab Adaptive \cr
#' \code{rk5} \tab \code{runge_kutta_dopri5} \tab Regular \cr
#' \code{rk5_a} \tab \code{runge_kutta_dopri5} \tab Adaptive \cr
#' \code{rk5_i} \tab \code{runge_kutta_dopri5} \tab Interpolating adaptive \cr
#' \code{rk78} \tab \code{runge_kutta_fehlberg78} \tab Regular \cr
#' \code{rk78_a} \tab \code{runge_kutta_fehlberg78} \tab Adaptive \cr
#' \code{abN} \tab \code{adams_bashforth} \tab Order N multistep \cr
#' \code{abmN} \tab \code{adams_bashforth_moulton} \tab Order N multistep \cr
#' \code{bs} \tab \code{bulirsch_stoer} \tab Adaptive \cr
#' \code{bsd} \tab \code{bulirsch_stoer_dense_out} \tab Interpolating adaptive}
#' 
#' These steppers are described at \href{http://headmyshoulder.github.io/odeint-v2/doc/boost_numeric_odeint/odeint_in_detail/steppers.html#boost_numeric_odeint.odeint_in_detail.steppers.stepper_overview}{here}.
#' 
#' @note The c++11 plugin is enabled.
#'  
#' @return The C++ code invisibly.
#' 
#' The following functions are generated:
#' 
#' \tabular{llll}{
#'  Function \tab Use \tab Arguments \tab Return \cr
#'  \code{name} \tab
#'    regular observer calls \tab
#'    \code{init, duration, step_size = 1.0, start = 0.0} \tab
#'    data frame\cr
#'  \code{name_adap} \tab
#'    adaptive observer calls \tab
#'    \code{init, duration, step_size = 1.0, start = 0.0} \tab
#'    data frame \cr
#'  \code{name_at} \tab
#'    specified observer calls \tab
#'    \code{init, times, step_size = 1.0} \tab
#'    data frame \cr
#'  \code{name_continue_at} \tab
#'    specified observer calls starting from previous final state \tab
#'    \code{times, step_size = 1.0} \tab
#'    data frame \cr
#'  \code{name_no_record} \tab
#'    no observer calls \tab
#'    \code{init, duration, step_size = 1.0, start = 0.0} \tab
#'    vector (final state)\cr
#'  \code{name_reset_observer} \tab
#'    clear observed record \tab void \tab void \cr
#'  \code{name_get_state} \tab
#'    get current state \tab void \tab vector \cr
#'  \code{name_set_state} \tab
#'    set current state \tab \code{new_state} \tab void \cr
#'  \code{name_get_output} \tab
#'    fetch observed record \tab void \tab data frame \cr
#'  \code{name_get_params} \tab
#'    get parameter values \tab void \tab a list \cr
#'  \code{name_set_params} \tab
#'    set parameter values \tab parameters \tab void}
#'    
#' Arguments are:
#' 
#' \tabular{ll}{
#' \code{init} \tab vector of initial conditions \cr
#' \code{duration} \tab end at start + duration \cr
#' \code{step_size} \tab the integration step size; variable for adaptive methods \cr
#' \code{start} \tab the starting time (always equal 0.0 for \code{name_at}) \cr
#' \code{time} \tab vector of times as which to call the observer \cr
#' \code{new_state} \tab vector of state values \cr}
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{set_optimization}}, \code{\link{integrate_sys}}
#' 
#' @examples
#' \dontrun{
#' # Logistic growth
#' compile_sys("logistic", "dxdt = x * (1 - x)")
#' plot(logistic(0.001, 15, 0.1), type = "l", lwd = 2, col = "steelblue")
#' Sys.sleep(0.5)
#' 
#' # Lotka-Volterra predator-prey equations
#' pars = c(alpha = 1, beta = 1, gamma = 1, delta = 1)
#' LV.sys = '
#'   dxdt[0] = alpha * x[0] - beta * x[0] * x[1];
#'   dxdt[1] = gamma * x[0] * x[1] - delta * x[1];
#' ' # LV.sys
#' compile_sys("preypred", LV.sys, pars, TRUE)
#' x = preypred(rep(2, 2), 100, 0.01)
#' plot(x[, 2:3], type = "l", lwd = 2,
#'      xlab = "Prey", ylab = "Predator",
#'      col = "steelblue")
#' Sys.sleep(0.5)
#' 
#' # Lorenz model from odeint examples
#' pars = c(sigma = 10, R = 28, b = 8 / 3)
#' Lorenz.sys = '
#'   dxdt[0] = sigma * (x[1] - x[0]);
#'   dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
#'   dxdt[2] = -b * x[2] + x[0] * x[1];
#' ' # Lorenz.sys
#' compile_sys("lorenz", Lorenz.sys, pars, TRUE)
#' x = lorenz(rep(1, 3), 100, 0.001)
#' plot(x[, c(2, 4)], type = 'l', col = "steelblue")
#' }
#' @export
compile_sys = function(name, sys,
                       pars = NULL,
                       const = FALSE,
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
  if (!is.null(pars))
  {
    pcode = handle_pars(pars, const)
    globals = paste0(pcode$globals, globals, collapse = ";\n")
    if (!is.null(pcode$getter))
      footers = paste0(pcode$getter, footers, collapse = "\n")
    if (!is.null(pcode$setter))
      footers = paste0(pcode$setter, footers, collapse = "\n")
  }
  sys = paste0(sys, collapse = "; ")
  stepper = make_stepper_type(method)
  stepper_constr = make_stepper_constr(method, atol, rtol)
  if (ceiling(sys_dim) < 1) sys_dim = get_sys_dim(sys)
  if (is.na(sys_dim))
  {
    sys = vectorize_1d_sys(sys)
    sys_dim = 1L
  }
  code = read_template("rcpp_template")
  code = gsub("__STEPPER_TYPE__", stepper, code)
  code = gsub("__STEPPER_CONSTRUCT__", stepper_constr, code)
  code = sub("__SYS_SIZE__", ceiling(sys_dim), code)
  code = sub("__GLOBALS__", globals, code)
  code = sub("__SYS__", sys, code)
  code = sub("__HEADERS__", headers, code)
  code = sub("__FOOTERS__", footers, code)
  code = gsub("__FUNCNAME__", name, code)
  code = paste0(code, collapse = "\n")
  if (compile)
    tryCatch(Rcpp::sourceCpp(code = code, ...),
             error = function(e) e)
  return(invisible(code))
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
