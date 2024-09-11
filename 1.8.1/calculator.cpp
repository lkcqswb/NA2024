#include <iostream>
#include <math.h>
double newton(double (*func)(double),double (*f1)(double),double x_0,int iter,double epsilon){
    double x=x_0;
    double u;
    for(int i=0;i<iter;i++){
        u=func(x);
        std::cout<<"iteration "<<i+1<< ' ' <<f1(x); 
        if(std::abs(u)<epsilon) break;
        x=x-u/f1(x);
        std::cout<< ' ' <<x<<std::endl;
    }
    return x;
}

double func(double x){
    return (4*pow(x,3)-2*x*x+3);
}
double f1(double x){
    return 12*x*x-4*x;
}


int main(){
    double x=newton(func,f1,1,10,1e-10);
    std::cout<<func(x);
}