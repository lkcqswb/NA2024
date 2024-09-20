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
    double a,b,epsilon,delta;
    int max_iteration;
public:
    bisection_method(const Function f,double a,double b,double epsilon,double delta,int max_iteration):
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


class Newton_method: public EquationSolver{
private:
    double x_0,epsilon,root;
    int max_iteration;
public:
    Newton_method(const Function f,double x_0,double epsilon,int max_iteration):
        EquationSolver(f),x_0(x_0),epsilon(epsilon),max_iteration(max_iteration)
    {}

    virtual double solve(){
        double diff,value;
        root=x_0;
        for (int i = 0; i < max_iteration; i++)
        {
            diff=F.differential(root);
            if(abs(diff)<1e-3) break;
            value=F.function_value(root);
            if(abs(value)<epsilon) break;
            root-=value/diff;
        }
        return root;    
    }

};

class Secant_method: public EquationSolver{
private:
    double x_0,x_1,epsilon,delta;
    int max_iteration;
public:
    Secant_method(const Function f,double x_0,double x_1,double epsilon,double delta,int max_iteration):
        EquationSolver(f),x_0(x_0),x_1(x_1),epsilon(epsilon),delta(delta),max_iteration(max_iteration)
    {}

    virtual double solve(){
        double diff,value_current,value_next,x_current=x_0,x_next=x_1;
        for (int i = 0; i < max_iteration; i++)
        {
            value_next=F.function_value(x_next);
            if(abs(value_next)<epsilon||abs(x_current-x_next)<delta) break;
            value_current=F.function_value(x_current);
            diff=(value_next-value_current)/(x_next-x_current);
            cout<<diff<<endl;
            if(abs(diff)<1e-6){exit(-1);}
            x_current=x_next;
            x_next-=value_next/diff;
            cout<<x_current<<" "<<x_next<<endl;
        }
        return x_next;    
    }

};