#include<iostream>
#include"EquationSolver.hpp"
#include<cmath>
double func(double x){
    return x-tan(x);
}

int main(){
    Function Function1(func);
    cout<<Newton_method(Function1,4.6,1e-3,10).solve();
}