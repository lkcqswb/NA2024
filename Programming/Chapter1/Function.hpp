class Function {
private:
    double (*function)(double);

public:
    Function(double (*f)(double)) : function(f) {}

    double function_value(double x) const {
        return function(x);
    }
};