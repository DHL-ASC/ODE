#include <nonlinfunc.h>
#include <ode.h>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};

class RCNetwork : public NonlinearFunction {
    double R = 1; // Ohm
    double C = 1e-3; // Farad

    size_t DimX() const override {
        return 2;
    }

    size_t DimF() const override {
        return 2;
    }
    
    void Evaluate (VectorView<double> x, VectorView<double> f) const override {
        f(0) = 1;
        f(1) = (std::cos(100 * M_PI * x(0)) - x(1)) / (R * C);
    }
    
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override {
        df = 0.0;
        df(0, 0) = 0;
        df(1, 1) = -1/(R*C);
        df(1, 0) = -1/(R*C) * 10 * M_PI * sin(10*M_PI*x(0));
    }
};



int main()
{
  double tend = 4*M_PI;
  int steps = 1000;
  Vector<> y(2); 
  y(0) = 1;
  y(1) = 0;
  auto rhs = make_shared<RCNetwork>();

  //SolveODE_EE(tend, steps, y, rhs,
  //            [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });

  //write steps from SolveODE_EE to file
  ofstream out("test_ode_ie.txt");
  SolveODE_IE(tend, steps, y, rhs,
              [&out](double t, VectorView<double> y) { out << t << "  " << y(0) << " " << y(1) << endl; }); 
  out.close();

}
