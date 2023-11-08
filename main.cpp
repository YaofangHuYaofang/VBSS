#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat computeSqrt(const arma::mat& inputMatrix) {
  arma::mat sqrtMatrix = arma::sqrt(inputMatrix); 
  
  return sqrtMatrix;
}

arma::vec computeEz(arma::vec Y, arma::vec E_z1, arma::vec E_z0) {
  int n = Y.size();
  
  arma::vec E_z = zeros(n);
  
  for (int i = 0; i < n; ++i) {
    if (Y[i] == 1)
      E_z(i, 0) = E_z1[i];  
    else if (Y[i] == 0)
      E_z(i, 0) = E_z0[i];  
  }
  
  return E_z;
}

arma::rowvec convertToRow(arma::vec p) {
  int n = p.size();
  
  arma::rowvec rowVector(n);
  
  for (int o = 0; o < n; ++o) {
    rowVector(o) = p(o);
  }
  
  return rowVector;
}

arma::mat combineVectorMatrix(const arma::vec& vec, const arma::mat& mat) {
  arma::mat combinedMatrix = arma::join_rows(mat, vec);  
  
  return combinedMatrix;
}

// [[Rcpp::export]]
Rcpp::List VBSS(int n, arma::mat A, arma::mat D,
                arma::vec key, arma::vec selected, arma::vec Y,
                double d, double tau, double sigma_a2, double sigma_b2,
                double epsilon_conv, int max_iter
){
  arma::mat inv_D = arma::inv(D); 
  arma::mat G = inv_D * A;
  
  // specify the k
  
  arma::mat w = computeSqrt(inv_D);
  mat B_star = diagmat(vec(n, fill::ones)) - d*w * A * w;
  vec theta_star = arma::inv(B_star) * Y;
  mat U = arma::inv(B_star) * trans(arma::inv(B_star));
  vec theta_0 = theta_star;
  
  mat P(n, 1, fill::zeros);
  P.fill(0.5);
  
  // set initial values for parameters sigma2, theta, alpha, a and b

  mat S_theta = E_sigma2*U;
  double E_theta2 = dot(E_theta, E_theta) + trace(S_theta);
  vec m = E_a*E_theta + E_b;
  vec t = zeros(n);
  for (int o = 0; o < n; ++o) {
    t(o) = R::pnorm(-m(o), 0.0, 1.0, true, false);
  }

  double E_logalpha = log(E_alpha);
  double g = trace(diagmat(vec(n, fill::ones)) / 2.0) + dot(E_theta - theta_0, (arma::inv(U) / 2.0) * (E_theta - theta_0));
  
  E_sigma2 = (tau + n/2)/(tau + g);
  
  vec E_z1 = zeros(n);
  for (int o = 0; o < n; ++o) {
    E_z1(o) = m[o] + R::dnorm(-m(o), 0.0, 1.0, false)/(1-t[o]);
  }
  
  vec part1 = zeros(n);
  
  for (int o = 0; o < n; ++o) {
    double numerator = std::sqrt(2 * M_PI) * m[o] * t[o] - std::exp(-std::pow(m[o], 2) / 2.0);
    double denominator = std::sqrt(2 * M_PI) *(t[o] + (1 - t[o]) * std::exp(E_logalpha));
    part1[o] = numerator / denominator;
  }
  
  vec part2 = zeros(n);
  
  for (int o = 0; o < n; ++o) {
    double numerator = std::exp(E_logalpha) *(std::sqrt(2 * M_PI) * m[o] * (1-t[o]) + std::exp(-std::pow(m[o], 2) / 2.0));
    double denominator = std::sqrt(2 * M_PI) * (t[o] + (1 - t[o]) * std::exp(E_logalpha));
    part2[o] = numerator / denominator;
  }
  
  vec E_z0 = part1 + part2;
  
  vec E_z = computeEz(Y, E_z1, E_z0);
  
  int converged = 0;
  
  for (int j = 1; j < max_iter; ++j){
    double S_a = 1/(E_theta2 + 1/sigma_a2);
    double tmp_E_a =  S_a * (dot(E_theta, E_z) - E_b * sum(E_theta));
    
    double S_b = sigma_b2 / (n * sigma_b2 + 1);
    double tmp_E_b = S_b * (sum(E_z) - tmp_E_a * sum(E_theta));

    mat tmp_S_theta = zeros(n, n);
    vec tmp_E_theta = zeros(n);
    
    for (int l = 0; l < n; ++l) {
      double S_theta_l = 1.0 / (S_b + E_sigma2 / U(l, l));
      double E_theta_l = S_theta_l * (tmp_E_b*E_z[l] - tmp_E_b*tmp_E_a + theta_0[l]*E_sigma2/U(l, l));
      
      tmp_S_theta(l, l) = S_theta_l;
      tmp_E_theta[l] = E_theta_l;
    }
    
    double tmp_E_theta2 = dot(tmp_E_theta, tmp_E_theta) + trace(tmp_S_theta);
    
    vec tmp_m = tmp_E_a*tmp_E_theta + tmp_E_b;
    vec tmp_t = zeros(n);
    for (int o = 0; o < n; ++o) {
      tmp_t(o) = R::pnorm(-tmp_m(o), 0.0, 1.0, true, false);
    }
    double tmp_g = trace(inv(U)*tmp_S_theta/2) + dot(tmp_E_theta - theta_0, (inv(U) / 2.0) * (tmp_E_theta - theta_0));
    
    vec tempvec1(n, 1); 
    for (int o = 0; o < n; ++o) {
      tempvec1[o] = exp(E_logalpha)*(1-Y[o])*(1-tmp_t[o])/(tmp_t[o]+exp(E_logalpha)*(1-tmp_t[o]));
    }
    
    double alpha_1 = sum(tempvec1) + 1.0;
    double beta_1 = sum(Y) + 1.0;
    
    double tmp_E_alpha = alpha_1 / (alpha_1 + beta_1);
    double tmp_E_logalpha = R::digamma(alpha_1) - R::digamma(alpha_1 + beta_1);
    
    double tmp_E_sigma2 = (tau + n / 2.0) / (tau + tmp_g);
    
    vec tmp_E_z1 = zeros(n);
    for (int o = 0; o < n; ++o) {
      tmp_E_z1(o) = tmp_m[o] + R::dnorm(-tmp_m(o), 0.0, 1.0, false)/(1-tmp_t[o]);
    }
    
    vec tmp_part1 = zeros(n);
    
    for (int o = 0; o < n; ++o) {
      double numerator = std::sqrt(2 * M_PI) * tmp_m[o] * tmp_t[o] - std::exp(-std::pow(tmp_m[o], 2) / 2.0);
      double denominator = std::sqrt(2 * M_PI) *(tmp_t[o] + (1 - tmp_t[o]) * std::exp(tmp_E_logalpha));
      tmp_part1[o] = numerator / denominator;
    }
    
    vec tmp_part2 = zeros(n);
    
    for (int o = 0; o < n; ++o) {
      double numerator = std::exp(tmp_E_logalpha) *(std::sqrt(2 * M_PI) * tmp_m[o] * (1-tmp_t[o]) + std::exp(-std::pow(tmp_m[o], 2) / 2.0));
      double denominator = std::sqrt(2 * M_PI) * (tmp_t[o] + (1 - tmp_t[o]) * std::exp(tmp_E_logalpha));
      tmp_part2[o] = numerator / denominator;
    }
    
    vec tmp_E_z = computeEz(Y, tmp_E_z1, E_z0);
    
    vec p = zeros(n);
    for (int o = 0; o < n; ++o) {
      p(o) = R::pnorm(tmp_m(o), 0.0, 1.0, true, false);
    }
    P = join_rows(P, p);
    
    if (mean(abs(P.col(j) - P.col(j-1))) < epsilon_conv) {
      Rcpp::Rcout << "Converged at iteration " << j << std::endl;
      converged = 1;
      break;
    }
    
    if (j == max_iter) { 
      Rcpp::warning("VB did not converge!\n"); 
    }
    
    if (j >= 1) {
      E_theta2 = tmp_E_theta2;
      E_theta = tmp_E_theta;
      E_z = tmp_E_z;
      E_b = tmp_E_b;
      E_a = tmp_E_a;
      E_sigma2 = tmp_E_sigma2;
      E_logalpha = tmp_E_logalpha;
    }
  }  
  
  return List::create(Named("probs") = P, Named("convergence") = converged);
}