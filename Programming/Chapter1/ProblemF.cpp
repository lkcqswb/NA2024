#include<iostream>
#include<cmath>
#include <functional>
#include"EquationSolver.hpp"


class nose_in_failure_solver
{
private:
    double l,h,D,beta1,beta1_pi;
public:
    nose_in_failure_solver(double l,double h,double D,double beta1):l(l),h(h),D(D),beta1(beta1){
    }
    
    void reset_param(double new_l,double new_h,double new_D,double new_beta1){
        l=new_l;
        h=new_h;
        D=new_D;
        beta1=new_beta1;
    }
    double get_param(string str){
        if(str=="l") return l;
        if(str=="h") return h;
        if(str=="D") return D;
        if(str=="beta1") return beta1;
        cout<<"invalid input"<<endl;
        exit(-1);
    }
};
nose_in_failure_solver solver1(89, 49, 55, 11.5);

double equation(double alpha){
        double alpha_pi = alpha * pi / 180,l=solver1.get_param("l"),h=solver1.get_param("h"),D=solver1.get_param("D"),beta1=solver1.get_param("beta1"),A,B,C,E;
        double beta1_pi=beta1*pi/180;
        A=l*sin(beta1_pi);
        B=l*cos(beta1_pi);
        C=(h+0.5*D)*sin(beta1_pi)-0.5*D*tan(beta1_pi);
        E=(h+0.5*D)*cos(beta1_pi)-0.5*D;
        return A * sin(alpha_pi) * cos(alpha_pi) + B * pow(sin(alpha_pi), 2) - C * cos(alpha_pi) - E * sin(alpha_pi);
    }





int main() {
    cout<<"Problem F"<<endl;
    Function func(&equation);
    cout<<"F-a"<<endl;
    cout <<"alpha:"<<Newton_method(func, 45, 1e-3, 200).solve() << endl;

    cout<<"F-b"<<endl;
    solver1.reset_param(89, 49, 30, 11.5);
    cout <<"alpha:"<< Newton_method(func, 33, 1e-3, 200).solve() << endl;


    cout<<"F-c"<<endl;
    //solver1.reset_param(89, 49, 55, 11.5);
    double inital_value[6]={10,70,80,87,150,140};
    for (int i = 0; i < 6; i+=2)
    {
        cout<<"using initial value:"<<"alpha_1="<<inital_value[i]<<" alpha_2="<<inital_value[i+1]<<endl;
        cout <<"alpha="<< Secant_method(func, inital_value[i], inital_value[i+1], 1e-3,1e-3,200).solve() << endl;
    }
    
    
    return 0;
}