#include<iostream>
#include<cmath>
int main(){
    double x=(6+sqrt(21))/5;
    std::cout<<x*pow(x-1,2)*pow(x-3,2)/120;
}