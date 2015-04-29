// [[Rcpp::export]]
std::vector<double> __FUNCNAME___get_params()
{
  return std::vector<double>(odeintr::pars.begin(), odeintr::pars.end());
}

