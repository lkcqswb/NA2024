#ifndef FUNCTION_H
#define FUNCTION_H
#define pi 3.1415926535
#define M_e 2.71828
class Function{
private:
    double (*function)(double);
public:
    Function(double (*func)(double)) : function(func) {}
    double function_value(double x) const {
        return function(x);
    }
};
#endif
