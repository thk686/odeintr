## Test environments
* local OS X install, R 3.3.2
* ubuntu 12.04 (on travis-ci), R 3.3.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* No reverse depends that I know of.

* Given that hundreds of CRAN packages generate C-call points using Rcpp, would it not make sense to coordinate with Rcpp folks on automating the registration of the C-functions prior to adding a CRAN policy check?

* I have attempted to work around problems related to testthat, compiling C++11 code and r-devel. I am not sure I have fixed it, but perhaps this will work better.

* I think the solution is just to remove the tests from the package using .Rbuildignore

* Fixed a regexp in .Rbuildignore