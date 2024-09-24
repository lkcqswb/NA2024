#include<iostream>
#include<cmath>
#include"EquationSolver.hpp"

using namespace std;
double f1(double x){
    return sin(x/2)-1;
}

double f2(double x){
    return pow(M_e,x)-tan(x);
}

double f3(double x){
    return pow(x,3)-12*pow(x,2)+3*x+1;
}

int main(){
    cout<<"Problem D"<<endl;
    Function F1(&f1),F2(&f2),F3(&f3);
    cout<<"Using the initial value that was provided"<<endl;
    cout<<"root of sin(x/2)-1 with initial value (0,pi/2): "<<Secant_method(F1,0,pi/2,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of e^x-tan(x) with initial value (1,1.4): "<<Secant_method(F2,1,1.4,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of x^3-12*x^2+3*x+1 with initial value (0,-0.5): "<<Secant_method(F3,0,-0.5,1e-4,1e-3,100).solve()<<endl;
    cout<<"Using the other initial value"<<endl;
    double init[2];
    init[0]=1;
    init[1]=pi/2;
    cout<<"root of sin(x/2)-1 with initial value ("<<init[0]<<','<<init[1]<<"): "<<Secant_method(F1,init[0],init[1],1e-4,1e-3,100).solve()<<endl;
    init[0]=1;
    init[1]=2;
    cout<<"root of e^x-tan(x) with initial value ("<<init[0]<<','<<init[1]<<"): "<<Secant_method(F2,init[0],init[1],1e-4,1e-3,100).solve()<<endl;
    init[0]=3;
    init[1]=-0.5;
    cout<<"root of x^3-12*x^2+3*x+1 with initial value ("<<init[0]<<','<<init[1]<<"): "<<Secant_method(F3,init[0],init[1],1e-4,1e-3,100).solve()<<endl;
}