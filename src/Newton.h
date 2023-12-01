#ifndef Newton_h
#define Newton_h

#include "nonlinfunc.h"
#include <math.h>
namespace ASC_ode
{

  void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<> res(func->DimF());
    Matrix<> fprime(func->DimF(), func->DimX());
    Matrix<> invfprime(func->DimF(), func->DimX());

    for (int i = 0; i < maxsteps; i++)
      {
        func->Evaluate(x, res);
        func->EvaluateDeriv(x, fprime);
        invfprime = fprime.Inverse();
        x -= invfprime*res;

        double err= sqrt(L2NormSquared(res));
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
