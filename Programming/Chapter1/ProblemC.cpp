#include<iostream>
#include"EquationSolver.hpp"
#include<cmath>
double func(double x){
    return x-tan(x);
}

int main(){
    cout<<"Problem C"<<endl;
    Function Function1(func);
    cout<<"root of x-tan(x) with initial value 4.6"<<Newton_method(Function1,4.6,1e-3,10).solve()<<endl;
}