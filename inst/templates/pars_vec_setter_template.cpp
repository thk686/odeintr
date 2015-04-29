// [[Rcpp::export]]
void __FUNCNAME___set_params(std::vector<double> p)
{
  if (p.size() != odeintr::pars.size())
    Rcpp::stop("Invalid parameter vector");
  std::copy(p.begin(), p.end(), odeintr::pars.begin());
}

