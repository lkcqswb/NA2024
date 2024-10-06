#include<iostream>
#include"EquationSolver.hpp"
#include<cmath>
double func(double x){
    return x-tan(x);
}

int main(){
    cout<<"Problem C"<<endl;
    Function Function1(func);
    cout<<"root of x-tan(x) with initial value 4.5: "<<Newton_method(Function1,4.5,1e-3,10).solve()<<endl;
    cout<<"root of x-tan(x) with initial value 7.7: "<<Newton_method(Function1,7.7,1e-3,10).solve()<<endl;
}