#include<iostream>
#include<cmath>
#include"EquationSolver.hpp"
#define pi 3.1415926535

class nose_in_failure_solver
{
private:
    double l,h,D,beta1,beta1_pi,A,B,C,E;
public:
    nose_in_failure_solver(double l,double h,double D,double beta1):l(l),h(h),D(D),beta1(beta1){
        beta1_pi=beta1*pi/180;
        A=l*sin(beta1_pi);
        B=l*cos(beta1_pi);
        C=(h+0.5*D)*sin(beta1_pi)-0.5*D*tan(beta1_pi);
        E=(h+0.5*D)*cos(beta1_pi)-0.5*D;
    }
    double equation(double alpha){
        double alpha_pi = alpha * pi / 180;
        return A * sin(alpha_pi) * cos(alpha_pi) + B * pow(sin(alpha_pi), 2) - C * cos(alpha_pi) - E * sin(alpha_pi);
    }
};





int main() {
    nose_in_failure_solver solver1(89, 49, 55, 11.5);

    cout << Newton_method(Function(&solver1.equation), 45, 1e-3, 200).solve() << endl;

    return 0;
}