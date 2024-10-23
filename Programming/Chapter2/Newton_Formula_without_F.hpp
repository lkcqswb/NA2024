#ifndef NEWTONFORMULAW_H
#define NEWTONFORMULAW_H
#include<iostream>
#include<vector>
using namespace std;
typedef unsigned long long size_t;
class  Newton_formula_WF{
protected:
    vector<vector<double>> interpolation_points;
    vector<double> coefficients;
public:
    Newton_formula_WF(vector<vector<double>> points):interpolation_points(points) {
        size_t i,j;
        for (i = 0; i < interpolation_points.size(); i++)
        {
            coefficients.push_back(interpolation_points[i][1]);
        }
        for (i = 1; i < interpolation_points.size(); i++)
        {
            for (j = interpolation_points.size()-1; j >=i; j--)
            {
                coefficients[j]=(coefficients[j]-coefficients[j-1])/(interpolation_points[j][0]-interpolation_points[j-i][0]);
            }
            
        }   
    }
    double get_interpolation_value(double x){
        double po=1,value=0;
        size_t i;
        for(i=0;i<coefficients.size();i++){
            value+=coefficients[i]*po;
            po*=x-interpolation_points[i][0];
        }
        return value;
    }
    
};

#endif