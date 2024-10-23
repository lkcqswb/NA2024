#include<iostream>
#include<vector>
#include<cmath>
using namespace std;
#define pi 3.1415926535
class Cubic_Bezier
{
private:
    double eps=1e-6;
    int m,dimension;
    vector<double> t;
    vector<vector<double>> points;
    vector<vector<vector<double>>> qs;
public:
    Cubic_Bezier(vector<vector<double>> interpolation_points,vector<double>  ts,vector<vector<vector<double>>> derivatives,int input_dimension){
        size_t i;
        int j;
        dimension=input_dimension;
        m=(int)(interpolation_points.size())-1;
        points=interpolation_points;
        t=ts;
        if((int)derivatives.size()<m+1){
            throw "Insufficient infomation";
        }
        for (i = 0; i < ts.size()-1; i++)
        {
            vector<vector<double>> q;
            q.push_back(points[i]);
            vector<double> new_point=points[i];
            for (j = 0; j < dimension; j++)
            {
                new_point[j]+=(derivatives[i][1][j])/3;
            }
            q.push_back(new_point);
            new_point.clear();
            new_point=points[i+1];
            for (j = 0; j < dimension; j++)
            {
                new_point[j]-=(derivatives[i+1][0][j])/3;
            }
            q.push_back(new_point);
            new_point.clear();
            q.push_back(points[i+1]);

            qs.push_back(q);
        }      
    }
    vector<double> get_value(double t_input){
        int location=0,i;
        for (location = 0; location <= m; location++)
        {
            if(t[location]>t_input) break;
        }
        if(location==0){
            throw "Out of range";
        }else if(location==m+1){
            if(t[m]+eps>t_input) location=m;
            else throw "Out of range";
        }
        double t_prime=(t_input-t[location-1])/(t[location]-t[location-1]);
        vector<double> result;
        for ( i = 0; i < dimension; i++)
        {
            result.push_back(qs[location-1][0][i]*pow(1-t_prime,3)+3*qs[location-1][1][i]*t_prime*pow(1-t_prime,2)+3*qs[location-1][2][i]*pow(t_prime,2)*(1-t_prime)+qs[location-1][3][i]*pow(t_prime,3));
        }

        return result;
    }
    
};




