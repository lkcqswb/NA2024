#include<iostream>
#include<cmath>
#include"EquationSolver.hpp"
double L=10,r=1,V=12.4;

double volumn(double h){
    return L*(0.5*pi*pow(r,2)-pow(r,2)*asin(h/r)-h*sqrt(pow(r,2)-pow(h,2)))-V;
}


int main(){
    Function Vol(&volumn);
    cout<<bisection_method(Vol,0,r,1e-2,1e-2,200).solve()<<endl;
    cout<<Newton_method(Vol,0.5*r,1e-2,200).solve()<<endl;
    cout<<Secant_method(Vol,0,r,1e-2,1e-2,200).solve()<<endl;
}