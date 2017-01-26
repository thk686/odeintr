// Copyright Timothy H. Keitt 2015

#ifndef __UTILS_H__
#define __UTILS_H__

static const int M = std::sqrt(N);

static int
bound(int x){
  int res = x % M;
  return res < 0 ? M + res : res;}

static double
d_left(const state_type& x, const int i, const int j){
  return x[i * M + bound(j - 1)] - x[i * M + j];}

static double
d_up(const state_type& x, const int i, const int j){
  return x[bound(i - 1) * M + j] - x[i * M + j];}

static double
d_right(const state_type& x, const int i, const int j){
  return x[i * M + bound(j + 1)] - x[i * M + j];}

static double
d_down(const state_type& x, const int i, const int j){
  return x[bound(i + 1) * M + j] - x[i * M + j];}

static double
laplace2D(const state_type& x, const int i, const int j){
  return d_left(x, i, j) + d_up(x, i, j) + d_right(x, i, j) + d_down(x, i, j);}

#endif // __UTILS_H__
