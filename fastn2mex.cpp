#include <limits>
#include <algorithm>
#include <memory>

#include "mex.h"
#include "mat.h" 

#define EPSILON   1e-13 // used in deciding the init interval
#define TOLERANCE 1e-14 // used in optimisation
#define MAXITER   1e3

// #define DEBUG

#define LOG_FLOAT(X) mexPrintf("[Line:%d], [%s:%f]\n", __LINE__, #X, X);
#define LOG_INT(X) mexPrintf("[Line:%d], [%s:%d]\n", __LINE__, #X, X);

#include "fzero.cpp"  

/*
  BLAS and LAPACK function
 */
extern "C" {
  void dpofa_(double*, int*, int*, int*);
  void dposl_(double*, int*, int*, double*);
  void dtrsm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
  void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
  void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
  void dsysv_(char*, int*, int*, double*, int*, int*, double*, int*, double*, int*, int*);
}

/* 
   Define the target function
*/

struct f_params
{
  int n; // number of terms.
  double* d; // d vector
  double* c2; // c^2 vector
};

double f_left(double x, const void *params) {
  const struct f_params *p = static_cast<const f_params*>(params);
  const double *d = p->d;
  const double *c2 = p->c2;

  double ret = 0;
  for (int i=0; i<p->n; i++) {
    ret += c2[i]*d[i]/pow(1+x*d[i],2);
  }
  return ret;
}

double f_right(double x, const void *params) {
  const struct f_params *p = static_cast<const f_params*>(params);
  const double *d = p->d;
  const double *c2 = p->c2;

  double ret = 0;
  for (int i=0; i<p->n; i++) {
    ret += c2[i]*d[i]/pow(x+d[i],2);
  }
  return ret;
}

enum WHICH_INTERVAL {INTERVAL_LEFT, INTERVAL_RIGHT};
double optim(double *c2, double *e, int n, double x_lo, double x_hi, WHICH_INTERVAL interval, int& info) {
  double r;
  struct f_params params;
  params.n = n;
  params.d = e;
  params.c2 = c2;
  void* pp = static_cast<void*>(&params);

  if (interval == INTERVAL_LEFT) {
    r = fzero(f_left, pp, x_lo, x_hi, TOLERANCE, MAXITER, info);
  }
  else {
    r = 1.0/fzero(f_right, pp, x_lo, x_hi, TOLERANCE, MAXITER, info);
  }

  return r;
}


void computeOptimX(int n, const double *mu, const double *S, const double &lambda,
		   const double &d, double *x, int& info) {
  // process sigma
  auto sigma = std::unique_ptr<double[]>(new double[n*n]);
  double prod = lambda*d;
  for (int i=0; i<n*n; i++) sigma[i] = S[i] - prod;
  for (int i=0; i<n*n; i+=n) {
    sigma[i] += lambda;
    i++;
  }
  // process x
  for (int i=0; i<n; i++) x[i] = mu[i];

  // compute
  char charU = 'U';
  int iONE = 1;
  auto ipiv = std::unique_ptr<int[]>(new int[n]);
  double optwlen;
  int lwork = -1;
  dsysv_(&charU, &n, &iONE, sigma.get(), &n, ipiv.get(), x, &n, &optwlen, &lwork, &info);
  lwork = static_cast<int>(optwlen);
  auto work = std::unique_ptr<double[]>(new double[lwork]);
  dsysv_(&charU, &n, &iONE, sigma.get(), &n, ipiv.get(), x, &n, work.get(), &lwork, &info);
  if (info!=0) {
    mexPrintf("Error in dsysv [line %d], [errcode: %d], [lam: %f]", __LINE__, info, lambda);
  }
}

double computeObj(const int& n, const double *mu, const double *S, const double*w) {
  double nomi = 0.0;
  double deno = 0.0;
  // compute the denominator
  double tmp = 0.0;
  for (int i=0; i<n; i++) {
    tmp = 0.0;
    for(int j=0; j<n; j++)
      tmp += S[i*n+j]*w[j];
    deno += tmp*w[i];
  }
  deno = std::sqrt(deno);
  // compute the nominator
  for (int i=0; i<n; i++) {
    nomi += mu[i]*w[i];
  }
  // return
  return nomi/deno;
}

/*
  Main routine
  - Input: prhs = [mu, S, d]
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // check flag: can be removed for speed
  if (nrhs != 3) {
    mexErrMsgTxt("should be 3 inputs");
  }

  // Sanity check: can be removed for speed
  int n = mxGetM(prhs[1]);
  if (std::max(mxGetM(prhs[0]), mxGetN(prhs[0]))!=n) {
    mexErrMsgTxt("dims of arguments 1 and 2 do not match");
  }

  /* 
     Do the work 
  */
  plhs[0]     = mxCreateDoubleMatrix(n,1,mxREAL);
  double *w   = mxGetPr(plhs[0]);
  
  double *mu  = mxGetPr(prhs[0]);
  double *S   = mxGetPr(prhs[1]);
  double d    = mxGetPr(prhs[2])[0];

  double lb_d = 1.0/n;

  // check special cases
  if (d < lb_d) {
    // not feasible, return NaN
    for (int i=0; i<n; i++) w[i] = std::numeric_limits<double>::quiet_NaN();
#ifdef DEBUG
    mexPrintf("NOT FEASIBLE\n");
#endif
    return;
  } else if (d == lb_d) {
    // only one solution
    for (int i=0; i<n; i++) w[i] = 1.0/n;
#ifdef DEBUG
    mexPrintf("UNIQUE SOLUTION\n");
#endif
    return;
  };

  // check if constraint is binding
  int info;

  auto L = std::unique_ptr<double[]>(new double[n*n]);
  for (int i=0; i<n*n; i++) L[i] = S[i];

  auto x = std::unique_ptr<double[]>(new double[n]);
  for (int i=0; i<n; i++) x[i] = mu[i];

  // x = S^{-1}mu
  // only the upper triangle of L is meaningful
  dpofa_(L.get(),&n,&n,&info);
  if (info!=0) {
    mexPrintf("Error in dpofa [line %d], [errcode: %d]", __LINE__, info);
  }
  dposl_(L.get(),&n,&n,x.get());

  // normalise x
  double sumx = 0.0;
  for (int i=0; i<n; i++) sumx += x[i];
  for (int i=0; i<n; i++) x[i] /= sumx;

  // check feasibility of x
  double norm2x = 0.0;
  for (int i=0; i<n; i++) norm2x += pow(x[i],2);

  // if x is feasible, return x
  if (sumx>0 && norm2x<=d) {
    for (int i=0; i<n; i++) w[i] = x[i];
#ifdef DEBUG
    mexPrintf("UNCONSTRAINED SOLUTION\n");
#endif    
    return;
  }

  // Check if maximum expected return is negative
  double a = 0.0;
  for (int i=0; i<n; i++) a += pow(mu[i],2);
  double b = 0.0;
  for (int i=0; i<n; i++) b += mu[i];
  double c = 1/lb_d;
  double maxmu = (b+std::sqrt((a*c-pow(b,2))*(c*d-1)))/c;
  if (maxmu < 0) {
    // TODO: return an error flag
  }

  // now do the optimisation
  double ONE = 1.0;
  //  double ZERO = 0.0;
  int iONE = 1;
  char charL = 'L';
  char charU = 'U';
  char charN = 'N';
  char charT = 'T';
  char charR = 'R';
  char charV = 'V';

  auto LAmdL = std::unique_ptr<double[]>(new double[n*n]); // L'^{-1}*(A-d)*L^{-1}
  // 1. let LAmdL=I-d
  std::fill_n(LAmdL.get(), n*n, -d);
  for (int i=0; i < n*n; i+=n) {
    LAmdL[i] += 1;
    i++;
  }

  // 2. let LAmdL=L'^{-1}*LAmdL
  dtrsm_(&charL, &charU, &charT, &charN, &n, &n, &ONE, L.get(), &n, LAmdL.get(), &n);

  // 3. let LAmdL=LAmdL*L^{-1}
  dtrsm_(&charR, &charU, &charN, &charN, &n, &n, &ONE, L.get(), &n, LAmdL.get(), &n);

  // 4. do eigenvalue decomposition
  auto E = std::unique_ptr<double[]>(new double[n]); // the eigenvalues
  int lwork = -1;
  double ows = 0;
  dsyev_(&charV, &charU, &n, LAmdL.get(), &n, E.get(), &ows, &lwork, &info); // obtain best size of working array
  lwork = static_cast<int>(ows);
  auto work = std::unique_ptr<double[]>(new double[lwork]);
  dsyev_(&charV, &charU, &n, LAmdL.get(), &n, E.get(), work.get(), &lwork, &info); // LAmdL now contains the eigenvectors
  // now LAmdL contains the orthogonal eigenvectors, that is, Q=LAmdL
  if (info!=0) {
    mexPrintf("Error in dsyev [line %d], [errcode: %d]\n", __LINE__, info);
  }

  // 5. compute c2
  auto c2 = std::unique_ptr<double[]>(new double[n]);
  for (int i=0; i<n; i++) c2[i] = mu[i];
  dtrsm_(&charL, &charU, &charT, &charN, &n, &iONE, &ONE, L.get(), &n, c2.get(), &n);
  auto tmpc2 = std::unique_ptr<double[]>(new double[n]);
  int idx = 0;
  for (int i=0; i<n; i++) {
    tmpc2[i] = 0.0;
    for (int j=0; j<n; j++) {
      tmpc2[i] += LAmdL[idx]*c2[j];
      idx ++;
    }
  }
  for (int i=0; i<n; i++) c2[i] = pow(tmpc2[i],2);

  // 6. optimise on the intervals
  double testL = 0.0;
  for (int i=0; i<n; i++) testL += x[i];
  testL = d*pow(testL,2);
  for (int i=0; i<n; i++) testL -= pow(x[i],2);

  double lamL, lamR, lambda, obj, objL, objR;
  double *sol;

  auto xL = std::unique_ptr<double[]>(new double[n]);
  auto xR = std::unique_ptr<double[]>(new double[n]);

  if (testL <= 0 ) {
    // need to search the left interval
    info = 0;
    lamL = optim(c2.get(), E.get(), n, 0, -1/E[0]-EPSILON, INTERVAL_LEFT, info);

    // test // right interval
    double testR = 0.0;
    for (int i=0; i<n; i++) testR += mu[i];
    testR = d/(n*d-1)*pow(testR,2);

    for (int i=0; i<n; i++) testR -= pow(mu[i],2);
    if (testR<0) {
      // search right interval
      info = 0;
      lamR = optim(c2.get(), E.get(), n, EPSILON, -E[0]-EPSILON, INTERVAL_RIGHT, info);
      mexPrintf("[line %d], [LamR:%f]\n", __LINE__, lamR);

      // compare the two interval
      computeOptimX(n, mu, S, lamR, d, xR.get(), info);
      computeOptimX(n, mu, S, lamL, d, xL.get(), info);

      objL = computeObj(n,mu,S,xL.get());
      objR = computeObj(n,mu,S,xR.get());

      if (objL > objR) {
	lambda = lamL;
	obj    = objL;
	sol    = xL.get();
      } else {
	lambda = lamR;
	obj    = objR;
	sol    = xR.get();
      }
    } else {
      // final lambda is from left interval
      lambda = lamL;
      computeOptimX(n, mu, S, lambda, d, xL.get(), info);
      obj    = computeObj(n, mu, S, xL.get());
      sol    = xL.get();
    }
  } else {
    // search the right interval
    // final lambda is from right interval
    info = 0;
    lambda = optim(c2.get(), E.get(), n, EPSILON, -E[0]-EPSILON, INTERVAL_RIGHT, info);
    mexPrintf("[line %d], [Lam:%f]\n", __LINE__, info);
    computeOptimX(n, mu, S, lambda, d, xR.get(), info);    
    obj    = computeObj(n, mu, S, xR.get());
    sol    = xR.get();
  }

  /* output */
  double sumsol = 0.0;
  for (int i=0; i<n; i++) sumsol += sol[i];
  for (int i=0; i<n; i++) w[i] = sol[i]/sumsol;

  return;
}
