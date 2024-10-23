#include<iostream>
#include<vector>
#include<cmath>
#include"Cubic_Bezier.hpp"
#include <fstream> 
using namespace std;
#define pi 3.1415926535
vector<double> norm(vector<double> x){
    double sum=sqrt(pow(x[0],2)+pow(x[1],2));
    return {x[0]/sum,x[1]/sum};
}

Cubic_Bezier draw(int m){
    int i;
    vector<double> t;
    vector<vector<double>> points;
    vector<vector<vector<double>>> derivative;
    double x_temp,y_temp,t_temp=0,dy,eps=1e-2;
    for ( i = 0; i < m+1; i++)
    {
        x_temp=sqrt(3)*cos(t_temp*2*pi);
        y_temp=2*(sqrt(sqrt(3)*abs(cos(2*pi*t_temp)))-sqrt(3)*sin(2*pi*t_temp))/3;
        if(t_temp>0.25+eps&&t_temp<0.75-eps){
            dy=2*(-2*pi*sqrt(3)*cos(2*pi*t_temp)+pi*sqrt(sqrt(3))*sin(2*pi*t_temp)/sqrt(abs(cos(2*pi*t_temp))))/3;
            derivative.push_back({{(-sqrt(3)*2*pi)*sin(t_temp*2*pi)/m,dy/m},{(-sqrt(3)*2*pi)*sin(t_temp*2*pi)/m,dy/m}});
        }else if(t_temp<0.25-eps||t_temp>0.75+eps){
            dy=2*(-2*pi*sqrt(3)*cos(2*pi*t_temp)-pi*sqrt(sqrt(3))*sin(2*pi*t_temp)/sqrt(abs(cos(2*pi*t_temp))))/3;
            derivative.push_back({{(-sqrt(3)*2*pi)*sin(t_temp*2*pi)/m,dy/m},{(-sqrt(3)*2*pi)*sin(t_temp*2*pi)/m,dy/m}});
        }else{
            dy=20;
            if(abs(t_temp-0.25)<eps) derivative.push_back({{(-sqrt(3)*2*pi)*sin(t_temp*2*pi)/m,dy/m},{(-sqrt(3)*2*pi)*sin(t_temp*2*pi)/m,-dy/m}});
            else derivative.push_back({{(-sqrt(3)*2*pi)*sin(t_temp*2*pi)/m,dy/m},{(-sqrt(3)*2*pi)*sin(t_temp*2*pi)/m,-dy/m}});
        }
        points.push_back({x_temp,y_temp});
        t.push_back(t_temp);
        
        t_temp+=1/(double)(m);
    }
    Cubic_Bezier Bezier(points,t,derivative,2);
    return Bezier;
}




int main(){
    int m,i,j;
    int m_to_paint[3]={10,40,60};

    for ( i = 0; i < 3; i++)
    {
        ofstream outfile("ProblemF.txt", ios::trunc);
        if (!outfile.is_open()) {
            cerr << "can't open ProblemF.txt" << endl;
            return 1;
        }
        m=m_to_paint[i];
        for ( j = 0; j < m+1; j++)
        {
            double t_temp=(double)(j)/m;
            outfile << sqrt(3)*cos(t_temp*2*pi) << " " << 2*(sqrt(sqrt(3)*abs(cos(2*pi*t_temp)))-sqrt(3)*sin(2*pi*t_temp))/3 << endl;
        }
        outfile << "#END#dots"<< endl;
        Cubic_Bezier C1=draw(m);
        for (j = 0; j < 101; j++) {
            vector<double> result=C1.get_value(0.01*j);
            outfile << result[0] << " " << result[1] << endl;
        }
        outfile << "#END# m=" <<m<< endl;
        string command = "python plot2.py ProblemF.txt";
        system(command.c_str());
    }
    

}