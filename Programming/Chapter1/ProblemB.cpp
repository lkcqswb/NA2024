#include<iostream>
#include"EquationSolver.hpp"
#include<cmath>
using namespace std;


double f1(double x){
    double epsilon=1e-6;
    if(x==0) x+=epsilon;
    return 1/x-tan(x);
}

double f2(double x){
    double epsilon=1e-6;
    if(x==0) x+=epsilon;
    return 1/x-pow(2,x);
}

double f3(double x){
    return pow(2,-x)+pow(M_e,x)+2*cos(x)-6;
}

double f4(double x){
    return (pow(x,3)+4*pow(x,2)+3*x+5)/(2*pow(x,3)-9*pow(x,2)+18*x-2);
}


int main(){
    Function Function1(&f1),Function2(&f2),Function3(&f3),Function4(&f4);
    cout<<"root of function1: "<<bisection_method(Function1,0,pi/2,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of function2: "<<bisection_method(Function2,0,1,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of function3: "<<bisection_method(Function3,1,3,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of function4: "<<bisection_method(Function4,0,4,1e-4,1e-3,100).solve()<<endl;
}