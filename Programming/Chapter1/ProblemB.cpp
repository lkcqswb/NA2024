#include<iostream>
#include"EquationSolver.hpp"
#include"Function.hpp"
#include<cmath>
#define pi 3.1415926535
using namespace std;


double f1(double x){
    double epsilon=1e-6;
    if(x==0) x+=epsilon;
    return 1/x-tan(x);
}

int main(){
    Function Function1=Function(&f1);
    cout<<bisection_method(Function1,0,pi/2,1e-3,1e-4,200).solve();
}