# odeintr
Timothy H. Keitt  
04/17/2015  

The odeintr is package for integrating differential equations in R. The integration engine is
the [Boost ODEINT package](http://www.odeint.com). Some features:

1. Simple specification of the ODE system
1. Intelligent defaults, easily overridden, used throughout
1. A wide range of intergration methods available for compiled system (see [stepper types](http://www.boost.org/doc/libs/1_58_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html#boost_numeric_odeint.odeint_in_detail.steppers.stepper_overview))
1. Fully automated compilation of ODE system specified in C++
1. Results returned as a simple data frame ready for analysis and plotting
1. Ability to specify a custom observer in R that can return aribtrary data (not yet for compiled code)
1. Three options for calling the observer: at regular intervals, after each update step or at specified times
1. Ability to alter system state and restart simulations where you left off
1. You can easily save and edit the generated C++ code (please respect copyright)

### Examples


```r
library(odeintr)
dxdt = function(x, t) x * (1 - x)
x = integrate_sys(dxdt, 0.001, 15, 0.01)
plot(x, type = "l", lwd = 3, col = "steelblue")
```

![](README_files/figure-html/unnamed-chunk-1-1.png) 

```r
compile_sys("logistic", "dxdt = x * (1 - x)")
plot(logistic(0.001, 15, 0.01), type = "l", lwd = "3", col = "steelblue")
points(logistic_at(0.001, sort(runif(10, 0, 15)), 0.01), col = "darkred")
```

![](README_files/figure-html/unnamed-chunk-1-2.png) 

```r
dxdt = function(x, t) c(x[1] - x[1] * x[2], x[1] * x[2] - x[2])
obs = function(x, t) c(Prey = x[1], Predator = x[2])
x = integrate_sys(dxdt, rep(2, 2), 20, 0.01, observer = obs)
plot(x[, c(2, 3)], type = "l", lwd = 2, col = "steelblue")
```

![](README_files/figure-html/unnamed-chunk-1-3.png) 

```r
# C++ code
Lorenz.sys = '
  dxdt[0] = 10.0 * (x[1] - x[0]);
  dxdt[1] = 28.0 * x[0] - x[1] - x[0] * x[2];
  dxdt[2] = -8.0 / 3.0 * x[2] + x[0] * x[1];
  ' # Lorenz.sys
compile_sys("lorenz", Lorenz.sys)
x = lorenz(rep(1, 3), 100, 0.001)
plot(x[, c(2, 4)], type = 'l', col = "steelblue")
```

![](README_files/figure-html/unnamed-chunk-1-4.png) 

### Performance

One of the main reasons the compiled code has the potential to be very fast is that ODEINT is a header-only library, so the entire integration path is exposed to the compiler. That means your system functions can be inlined with the integration code, loops unrolled, etc. It will help if you enable optimziation in your compiler. Use "-O3" with gcc. See the R documentation on the user Makevars file. (Odeintr now provides a convenient function to set the compiler
optimization level.)

### To Do

1. ~~Add additional integration methods from odeint~~
1. Extend customized observer to compiled code
1. ~~Allow user to set error tolerances in compiled code~~
1. Allow user to set error tolerances for system defined in R
1. Expose implicit solver methods

Pull requests are welcome.

### Installation

```
devtools::install_github("thk686/odeintr")
```
