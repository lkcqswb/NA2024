#ifndef NEWTONFORMULA_H
#define NEWTONFORMULA_H
#include<iostream>
#include"Function.hpp"
#include<vector>
using namespace std;

class  Newton_formula{
protected:
    const Function & F;
    vector<double> interpolation_points;
    vector<double> coefficients;
public:
    Newton_formula(const Function& function,vector<double> points):F(function),interpolation_points(points) {
        size_t i,j;
        for (i = 0; i < interpolation_points.size(); i++)
        {
            coefficients.push_back(F.function_value(interpolation_points[i]));
        }
        for (i = 1; i < interpolation_points.size(); i++)
        {
            for (j = interpolation_points.size()-1; j >=i; j--)
            {
                if((interpolation_points[j]==interpolation_points[j-i])) throw "Input error";
                coefficients[j]=(coefficients[j]-coefficients[j-1])/(interpolation_points[j]-interpolation_points[j-i]);
            }
            
        }   
    }
    double get_interpolation_value(double x){
        double po=1,value=0;
        size_t i;
        for(i=0;i<coefficients.size();i++){
            value+=coefficients[i]*po;
            po*=x-interpolation_points[i];
        }
        return value;
    }
    
};

#endif