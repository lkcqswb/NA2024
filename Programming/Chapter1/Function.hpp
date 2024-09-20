#include<iostream>
class Function{
private:
    double (*function)(double);
    double eps=1e-6;
public:
    Function(double (*func)(double)) : function(func) {}
    double function_value(double x) const {
        return function(x);
    }
    double differential(double x) const{
        return (function(x+eps)-function(x))/eps;
    }

};