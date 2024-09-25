#include<iostream>
#include"EquationSolver.hpp"
#include<cmath>
using namespace std;


double f1(double x){
    double epsilon=1e-6;
    if(x==0) x+=epsilon;
    if(x==pi/2) x-=epsilon;
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
    cout<<"Problem B"<<endl;
    Function Function1(f1),Function2(f2),Function3(f3),Function4(f4);
    cout<<"root of 1/x-tan(x) on [0,pi/2]: "<<bisection_method(Function1,0,pi/2,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of 1/x-2^x on [0,1]: "<<bisection_method(Function2,0,1,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of 2^-x+e^x+2*cos(x)-6 on [1,3]: "<<bisection_method(Function3,1,3,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of (x^3+4*x^2+3*x+5)/(2*x^3-9*x^2+18*x-2):"<<endl;
    cout<<"If direct use of the bisection method is employed, the answer is: "<<bisection_method(Function4,0,4,1e-4,1e-3,100).solve()<<endl;
    cout<<"However, it is noted that the denominator has a zero point within the interval, "<<endl;
    cout<<"and the derivative of the denominator is greater than 0, "<<endl;
    cout<<"with the denominator having different signs on either side of the zero point. "<<endl;
    cout<<"Meanwhile, the numerator is always greater than 0 within the interval, "<<endl;
    cout<<"so the solution obtained is actually the root of the denominator."<<endl;
    cout<<"The equation itself has no roots."<<endl;
}