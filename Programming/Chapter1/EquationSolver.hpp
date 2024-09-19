#include<iostream>
#include"Function.hpp"
using namespace std;
class  EquationSolver{
protected:
    const Function & F;
public:
    EquationSolver(const Function& function):F(function) {}
    virtual double solve()=0;
};

class bisection_method: public EquationSolver{
private:
    double a,b,epsilon,delta,max_iteration;
public:
    bisection_method(const Function f,double a,double b,double epsilon,double delta,double max_iteration):
        EquationSolver(f),a(a),b(b),epsilon(epsilon),delta(delta),max_iteration(max_iteration)
    {}

    virtual double solve(){
        double len=b-a,f_a,f_c;
        if(F.function_value(a)*F.function_value(b)>0){
            cout<<"precondition"<<endl;
            exit(-1);
        }
        for (int i = 0; i < max_iteration; i++)
        {
            len/=2;
            if(abs(f_c)<epsilon||len<delta) break;
            f_c=F.function_value(a+len);
            f_a=F.function_value(a);
            if(abs(f_c)<epsilon) break;
            a+=(f_a*f_c>0)*len;
        }
        return a+len;    
    }

};