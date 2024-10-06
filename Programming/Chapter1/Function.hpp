#define pi 3.1415926535
#define M_e 2.71828
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
        return (function(x+eps)-function(x-eps))/(2*eps);
    }

};