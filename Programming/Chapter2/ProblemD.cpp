#include<iostream>
#include"Hermite.hpp"
#include<cmath>
using namespace std;

int main(){
    vector<vector<double>> dots={{0,0,1},{3,225,1},{5,383,1},{8,623,1},{13,993,1}};
    vector<vector<double>> derivative={{75},{77},{80},{74},{72}};
    Hermite H1(dots,derivative);
    cout<<"Predict the position of the car and its speed for t=10 seconds"<<endl;
    cout<<"position: "<<H1.get_interpolation_value(10)<<" feets. "<<" speed: "<<H1.get_derivative_value(10)<<" feets per second."<<endl;
    
    double t=8,lr=0.1,eps=1e-6;//gradient descent
    int max_iter=100,iter=0;
    while(H1.get_derivative_value(t)<=81){
        t+=lr*(H1.get_derivative_value(t+eps)-H1.get_derivative_value(t-eps))/(2*eps);
        iter+=1;
        if(t>=13||t<=0||iter>max_iter) break;
    }
    cout<<"This car exceeded the speed limit."<<endl;
    cout<<"t="<<t<<"s "<<" speed = "<<H1.get_derivative_value(t)<<"feets per second"<<endl;;
}