#include<iostream>
#include<cmath>
#include <iomanip>
#include"EquationSolver.hpp"
double L=10,r=1,V=12.4;

double volumn(double h){
    return L*(0.5*pi*pow(r,2)-pow(r,2)*asin(h/r)-h*sqrt(pow(r,2)-pow(h,2)))-V;
}


int main(){
    Function Vol(&volumn);
    cout<<"Problem E"<<endl;
    cout<<"using bisection method: "<<fixed << setprecision(2) <<bisection_method(Vol,0,r,1e-3,1e-4,500).solve()<<endl;
    cout<<"using Newton method: "<<fixed << setprecision(2)<<Newton_method(Vol,0.5*r,1e-3,400).solve()<<endl;
    cout<<"using secant method: "<<fixed << setprecision(2)<<Secant_method(Vol,0,r,1e-3,1e-3,400).solve()<<endl;
}