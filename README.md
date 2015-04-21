# odeintr
Timothy H. Keitt  
04/17/2015  

The odeintr is package for integrating differential equations in R. The integration engine is
the [Boost ODEINT package](http://www.odeint.com). Some features:

1. You can specify your system in C++ and odeintr will compile and link a custom integrator
1. A function that will write a user Makevars file with specified compiler
optimization level
1. 1-dimensional systems can omit C++ array indexing
2. If you define your system in R code, you can provide a custom observer function that can record arbitrary information after each integration step

### Examples


```r
library(odeintr)
dxdt = function(x, t) x * (1 - x)
y = integrate_sys(dxdt, 0.001, 15)
plot(y, type = "l", lwd = 3, col = "steelblue")
```

![](README_files/figure-html/unnamed-chunk-1-1.png) 

```r
compile_sys("logistic", "dxdt = x * (1 - x)")
plot(logistic(0.001, 15), type = "l", lwd = "3", col = "steelblue")
```

![](README_files/figure-html/unnamed-chunk-1-2.png) 

```r
dxdt = function(x, t) c(x[1] - x[1] * x[2], x[1] * x[2] - x[2])
obs = function(x, t) c(Prey = x[1], Predator = x[2])
y = integrate_sys(dxdt, rep(2, 2), 20, observer = obs)
plot(y[, c(2, 3)], type = "l", lwd = 2, col = "steelblue")
```

![](README_files/figure-html/unnamed-chunk-1-3.png) 

```r
# C++ code
Lorenz.globals = '
  const double sigma_ = 10.0;
  const double R_ = 28.0;
  const double b_ = 8.0 / 3.0;
' # Lorenz.globals
Lorenz.sys = '
  dxdt[0] = sigma_ * (x[1] - x[0]);
  dxdt[1] = R_ * x[0] - x[1] - x[0] * x[2];
  dxdt[2] = -b_ * x[2] + x[0] * x[1];
  ' # Lorenz.sys
compile_sys("lorenz", Lorenz.sys, Lorenz.globals)
x = lorenz(rep(1, 3), 100)
plot(x[, c(2, 4)], type = 'l', col = "steelblue")
```

![](README_files/figure-html/unnamed-chunk-1-4.png) 

### Performance

One of the main reasons the compiled code has the potential to be very fast is that ODEINT is a header-only library, so the entire integration path is exposed to the compiler. That means your system functions can be inlined with the integration code, loops unrolled, etc. It will help if you enable optimziation in your compiler. Use "-O3" with gcc. See the R documentation on the user Makevars file. (Odeintr now provides a convenient function to set the compiler
optimization level.)

### To Do

1. Add additional integration methods from odeint
1. Extend customized observer to compiled code
1. Allow user to set tolerances

Pull requests are welcome.

### Installation

```
devtools::install_github("thk686/odeintr")
```
