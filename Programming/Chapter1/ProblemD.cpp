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
    Function F1(&f1),F2(&f2),F3(&f3);
    cout<<"root of Function1: "<<Secant_method(F1,0,pi/2,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of Function2: "<<Secant_method(F2,0,1.4,1e-4,1e-3,100).solve()<<endl;
    cout<<"root of Function3: "<<Secant_method(F3,0,-0.5,1e-4,1e-3,100).solve()<<endl;
}