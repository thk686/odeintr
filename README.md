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
1. You can easily save and edit the generated C++ code

### Examples


```r
library(odeintr)
dxdt = function(x, t) x * (1 - x)
system.time({x = integrate_sys(dxdt, 0.001, 15, 0.01)})
```

```
##    user  system elapsed 
##   0.118   0.005   0.123
```

```r
plot(x, type = "l", lwd = 3, col = "steelblue", main = "Logistic Growth")
```

![](README_files/figure-html/unnamed-chunk-1-1.png) 

```r
compile_sys("logistic", "dxdt = x * (1 - x)")
system.time({x = logistic(0.001, 15, 0.01)})
```

```
##    user  system elapsed 
##   0.000   0.000   0.001
```

```r
plot(x, type = "l", lwd = "3", col = "steelblue", main = "Logistic Growth")
points(logistic_at(0.001, sort(runif(10, 0, 15)), 0.01), col = "darkred")
```

![](README_files/figure-html/unnamed-chunk-1-2.png) 


```r
dxdt = function(x, t) c(x[1] - x[1] * x[2], x[1] * x[2] - x[2])
obs = function(x, t) c(Prey = x[1], Predator = x[2], Ratio = x[1] / x[2])
system.time({x = integrate_sys(dxdt, rep(2, 2), 20, 0.01, observer = obs)})
```

```
##    user  system elapsed 
##   0.215   0.017   0.233
```

```r
plot(x[, c(2, 3)], type = "l", lwd = 2, col = "steelblue", main = "Lotka-Volterra Phase Plot")
```

![](README_files/figure-html/unnamed-chunk-2-1.png) 

```r
plot(x[, c(1, 4)], type = "l", lwd = 2, col = "steelblue", main = "Prey-Predator Ratio")
```

![](README_files/figure-html/unnamed-chunk-2-2.png) 


```r
# C++ code
Lorenz.sys = '
  dxdt[0] = 10.0 * (x[1] - x[0]);
  dxdt[1] = 28.0 * x[0] - x[1] - x[0] * x[2];
  dxdt[2] = -8.0 / 3.0 * x[2] + x[0] * x[1];
  ' # Lorenz.sys
compile_sys("lorenz", Lorenz.sys)
system.time({x = lorenz(rep(1, 3), 100, 0.001)})
```

```
##    user  system elapsed 
##   0.012   0.000   0.012
```

```r
plot(x[, c(2, 4)], type = 'l', col = "steelblue", main = "Lorenz Attractor")
```

![](README_files/figure-html/unnamed-chunk-3-1.png) 


```r
VdP.sys = '
dxdt[0] = x[1];
dxdt[1] = 2.0 * (1 - x[0] * x[0]) * x[1] - x[0];
' # Vdp.sys
compile_sys("vanderpol", VdP.sys, method = "bsd") # Bulirsch-Stoer
system.time({x = vanderpol(rep(1e-4, 2), 100, 0.01)})
```

```
##    user  system elapsed 
##   0.003   0.000   0.002
```

```r
par(mfrow = c(2, 2), mar = rep(0.5, 4), oma = rep(5, 4), xpd = NA)
make.plot = function(xy, xlab = NA, ylab = NA)
  plot(xy, col = "steelblue", lwd = 2, type = "l",
       axes = FALSE, xlab = xlab, ylab = ylab)
plot.new()
make.plot(x[, c(3, 1)]); axis(3); axis(4)
make.plot(x[, c(1, 2)], "Time", "X1"); axis(1); axis(2)
make.plot(x[, c(3, 2)], "X2"); axis(1); axis(4)
title(main = "Van der Pol Oscillator", outer = TRUE)
```

![](README_files/figure-html/unnamed-chunk-4-1.png) 


```r
# Use a dynamic parameter
VdP.sys = '
dxdt[0] = x[1];
dxdt[1] = mu * (1 - x[0] * x[0]) * x[1] - x[0];
' # Vdp.sys
compile_sys("vpol2", VdP.sys, "mu", method = "bsd")
par(mfrow = c(2, 2), mar = rep(1, 4), oma = rep(3, 4), xpd = NA)
for (mu in seq(0.5, 2, len = 4))
{
  vpol2_set_params(mu = mu)
  x = vpol2(rep(1e-4, 2), 100, 0.01)
  make.plot(x[, 2:3]); box()
  title(paste("mu =", round(mu, 2)))
}
title("Van der Pol Oscillator Parameter Sweep", outer = TRUE)
title(xlab = "X1", ylab = "X2", line = 0, outer = TRUE)
```

![](README_files/figure-html/unnamed-chunk-5-1.png) 

### Performance

Because ODEINT is a header-only library, the entire integration path is exposed to the compiler. That means your system functions can be inlined with the integration code, loops unrolled, etc. It will help if you enable optimziation in your compiler. Use "-O3" with gcc. See the R documentation on the user Makevars file. (Odeintr now provides a convenient function to set the compiler
optimization level.)

The Lorenz  and Van der Pol examples above show about 10 million observer calls per second.

### To Do

1. ~~Add additional integration methods from odeint~~
1. ~~Extend customized observer to compiled code~~
1. ~~Allow user to set error tolerances in compiled code~~
1. Allow user to set error tolerances for system defined in R
1. Expose implicit solver methods
1. ~~Convenient dynamic parameter settings~~

Pull requests are welcome.

### Installation

```
devtools::install_github("thk686/odeintr")
```
