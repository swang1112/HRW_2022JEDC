#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// creat design matrix for VAR(p)
// [[Rcpp::export]]

arma::mat Ylag(arma::mat y, int p, int c ){
  
  arma::mat YLag = arma::zeros(y.n_rows, y.n_cols * p);
  
  for(int i = 0; i < p; ++i){
    YLag.submat(i, i * y.n_cols, y.n_rows-1, i * y.n_cols + (y.n_cols-1)) = y.submat(0, 0, y.n_rows - i - 1, y.n_cols-1 );
  }
  
  arma::mat YLagOut = YLag.submat(p-1,0, YLag.n_rows -2, YLag.n_cols - 1);
  
  if (c == 0)
  {
    return YLagOut;
  } else if (c == 1)
  {
    arma::vec c1(YLagOut.n_rows, arma::fill::ones);
    arma::mat YLagOut1 = arma::join_horiz(c1, YLagOut);
    return YLagOut1;
  } else if (c == 2)
  {
    arma::vec c1(YLagOut.n_rows, arma::fill::ones);
    arma::vec c2 = arma::linspace(1, YLagOut.n_rows, YLagOut.n_rows);
    arma::mat YLagOut2 = arma::join_horiz(c1, c2, YLagOut);
    return YLagOut2;
  } else if (c == 3)
  {
    arma::vec c1(YLagOut.n_rows, arma::fill::ones);
    arma::vec c2 = arma::linspace(1, YLagOut.n_rows, YLagOut.n_rows);
    arma::vec c3 = arma::pow(c2, 2);
    arma::mat YLagOut3 = arma::join_horiz(c1, c2, c3, YLagOut);
    return YLagOut3;
  } else
  {
    arma::mat place_holder = arma::zeros<arma::mat>(1,1);
    return place_holder;
  }
  
}

// LS estimation of unrestricted VAR model 
// y is vectorized 
// [[Rcpp::export]]

arma::mat var_est1(arma::mat y, arma::mat Z, int& K){
  
  arma::mat ZZ = arma::inv(Z*arma::trans(Z)) * Z; 
  arma::mat I  = arma::eye(K, K);
  arma::mat Out = arma::kron(I, ZZ) * y;
  return Out;
  
}

// matrix exponentials
arma::mat matexp(arma::mat X, int n)
{
  if (n == 0)
  {
    return arma::eye(X.n_cols, X.n_rows);
  } else if (n == 1)
  {
    return X;
  } else 
  {
    return X * matexp(X, n-1);
  }
}


// IRFs
// [[Rcpp::export]]
List IRF_fast(arma::mat &A_hat, arma::mat &B_hat, int &horizon)
{
  int K = A_hat.n_rows;
  int p = A_hat.n_cols / K;
  List Out(horizon);
  Out[0] = B_hat;
  if (p == 1)
  {
    for(int i = 1; i < horizon; i++)
    {
      Out[i] = matexp(A_hat, i) * B_hat;
    }
    return Out;
  } else
  {
    arma::mat Mm(K * p, K * p, arma::fill::zeros);
    Mm.submat(0, 0, (K - 1), (K * p - 1)) = A_hat;
    Mm.submat(K, 0, (K * p - 1), ((p - 1) * K - 1)) = arma::eye(K * (p - 1), K * (p - 1));
    arma::mat Mm1(K * p, K * p, arma::fill::eye);
    for (int i = 0; i < (horizon - 1); i++)
    {
      Mm1 = Mm1 * Mm;
      Out[i + 1] = Mm1.submat(0, 0, (K - 1), (K - 1)) * B_hat;
    }
    return Out;
  }
  
}


// draw from Inverse Wishart distribution
mat iwpq(int v, mat ixpx)
{
  int k = ixpx.n_rows;
  mat z = zeros(v,k);
  mat cixpx = chol(ixpx).t();
  for(int i=0; i<v; i++)
  {
   z.row(i) = trans(cixpx*randn(k,1));
  }
  return inv(z.t()*z);
}

  



// normal-inverse-wishart conjugate for a K-dimensional VAR(p) with an intercept
// [[Rcpp::export]]
List Gibbs_VAR_noninfo(arma::mat y_vec, arma::mat b_vec, arma::mat X, arma::mat Covmat,
                  int MC, int p, int K, int Tob)
{
  arma::mat Sigma = Covmat;
  List Out(MC);
  
  int pass;
  mat V;
  mat b;
  mat b_mat;
  mat A_L;
  mat A;
  cx_vec eigval;
  vec eigval_r;
  mat resid_mat;
  mat resid;
  mat scale;
  
  for (int m = 0; m < MC+1000; m++)
  {
    V = arma::inv_sympd(arma::kron(arma::inv_sympd(Sigma), X.t()*X));
    // draw AR coefficients
    pass = 0;
    b = b_vec + (arma::randn(1,b_vec.n_elem) * arma::chol(V)).t();
    while (pass == 0)
    {
      // turn b into a companion form and check stablility
      b_mat = reshape(b, K*p+1, K).t();
      // Rcout << b_mat.cols(1, b_mat.n_cols-1) << std::endl;
      A_L = join_horiz(eye<mat>(K*p-K, K*p-K), zeros<mat>(K*p-K, K));
      A   = join_vert(b_mat.cols(1, b_mat.n_cols-1), A_L);
      eigval = eig_gen( A ); 
      if (max(abs(eigval))<1)
      {
        pass = 1;
      } else
      {
        b = b_vec + (arma::randn(1,b_vec.n_elem) * arma::chol(V)).t();
        pass = 0;
      }
    }
    
    
    // draw Sigma
    pass = 0;
    resid = y_vec - kron(eye(K, K), X) * b;
    resid_mat = reshape(resid, X.n_rows, K);
    
    scale = resid_mat.t()*resid_mat + eye(K,K)*pow(Tob, -3);
    
    // Rcout << scale << std::endl;
    
    Sigma = iwpq(Tob, inv(scale));
    while (pass == 0)
    {
      // Rcout << Sigma << std::endl;
      eigval_r = eig_sym( Sigma ); 
      if (all(eigval_r > 0))
      {
        pass = 1;
      } else
      {
        Sigma = iwpq(Tob, inv(scale));
        pass = 0;
      }
    }
    
    
    if (m > 999){
      List Subout;
      Subout["b_vec"] = b;
      Subout["Sigma"] = Sigma;
      
      Out[m-1000] = Subout;    
    }
    
    
    //Out[m] = join_horiz(b_mat, Sigma);
  }
  
  
  return Out;
}
  

  
// point estimates proposed by Inoue and Kilian (2021)
// IRFs has dimension N times n_irf, where N = # draws
// loss_fun = 1 for quadratic and 2 for absolute
// [[Rcpp::export]]

vec IK_point(vec loss_val, mat IRFs, int loss_fun){
  int N = IRFs.n_rows;  
  double foo;
  if (loss_fun == 1)
  {
    for (int i=0; i<N; i++)
    {
      foo = 0.00;
      for (int m = 0; m < N; m++)
      {
        foo = foo + sum(pow(IRFs.row(i) - IRFs.row(m), 2));
      }  
      loss_val(i) = foo/N;
    }
    
  } else if (loss_fun == 2)
  {
    for (int i=0; i<N; i++)
    {
      foo = 0.00;
      for (int m = 0; m < N; m++)
      {
        foo = foo + sum(abs(IRFs.row(i) - IRFs.row(m)));
      }  
      loss_val(i) = foo/N;
    }
    
  } else
  {
  }
  
 return loss_val;
}