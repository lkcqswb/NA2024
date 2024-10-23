#ifndef HERMITE_H
#define HERMITE_H
#include<iostream>
#include"Function.hpp"
#include<vector>
using namespace std;
typedef unsigned long long size_t;

double factorial(int n){
    double sum=1;
    int i;
    if(n==0) return 1;
    for ( i = 1; i < n+1; i++)
    {
        sum*=i;
    }
    return sum;
}

class  Hermite{
protected:
    vector<vector<double>> interpolation_points;
    vector<vector<double>> derivative;
    vector<double> interpolation;
    vector<double> coefficients;
    void update(vector<double>* value,double x,int iter){
        int i;
        for (i = 0; i < (int)(coefficients.size()); i++)
        {
            if(iter==-1) value->push_back(1);
            else if(i!=iter) value->at(i)*=x-interpolation[iter];
        }
        
    }
public:
    Hermite(vector<vector<double>> points,vector<vector<double>> deriv):interpolation_points(points),derivative(deriv) {
        size_t i,j;
        for (i = 0; i < interpolation_points.size(); i++)
        {
           for(j=0;j<interpolation_points[i][2]+1;j++){
                coefficients.push_back(interpolation_points[i][1]);
                interpolation.push_back(interpolation_points[i][0]);
           }
        }
        for (i = 1; i < coefficients.size(); i++)
        {
            //output_coefficient();
            int pointer=interpolation_points.size()-1,value=interpolation_points[pointer][2];
            for (j = coefficients.size()-1; j >=i; j--)
            {
                if(interpolation[j]==interpolation[j-i]){
                    if(derivative[pointer].size()<i) throw "Insufficient differential information.";
                    coefficients[j]=derivative[pointer][i-1]/factorial(i);
                }
                else coefficients[j]=(coefficients[j]-coefficients[j-1])/(interpolation[j]-interpolation[j-i]);
                value-=1;
                if(value<0){
                    pointer-=1;
                    value=interpolation_points[pointer][2];
                }
            }   
        }   
    }
    double get_interpolation_value(double x){
        double po=1,value=0;
        size_t i;
        for(i=0;i<coefficients.size();i++){
            value+=coefficients[i]*po;
            po*=x-interpolation[i];
        }
        return value;
    }
    double get_derivative_value(double x){
        double result=0;
        vector<double> value; 
        size_t i,j;
        update(&value,x,-1);
        for(i=1;i<coefficients.size();i++){
            update(&value,x,(int)(i)-1);
            for ( j = 0; j < i; j++)
            {
                result+=coefficients[i]*value[j];
            }
        }
        return result;
    }
    void output_coefficient(){
        size_t i;
       for (i = 0; i < coefficients.size(); i++)
       {
         cout<<coefficients[i]<<" ";
       }
       cout<<endl; 
    }
    
};

#endif