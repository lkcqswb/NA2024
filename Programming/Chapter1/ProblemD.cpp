#include<iostream>
#include<cmath>
#include"EquationSolver.hpp"
#define pi 3.1415926535
using namespace std;
double f1(double x){
    return sin(x/2)-1;
}

int main(){
    Function F1=Function(&f1);
    cout<<Secant_method(F1,0,pi/2,1e-4,1e-3,100).solve();
}