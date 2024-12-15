#pragma once
#include<vector>
#include<iostream>
#include"../../include/Eigen/Dense"
#define pi 3.1415926535
using namespace std;
#include <cstdlib>
#include <ctime>
#include <cmath>

double sign(double x){
    if(x>=0) return 1;
    return -1;
}

void process_value_ud(vector<vector<double>> &dots, vector<double> centre,double r){
    static double angle=0;
    for (size_t i = 0; i < dots.size(); i++)
    {
        if(dots[i].size()<3){
            cout<<"invalid input"<<endl;
            throw "invalid input";
        }
    double x=dots[i][0]-centre[0],y=dots[i][1]-centre[1],z=dots[i][2]-centre[2];
    double r2=sqrt(x*x+y*y);
    if(r2==0){
        dots[i]= {asin(z/r),angle};
        continue;
    }
    double angle2=y>0 ? acos(x/r2):2*pi-acos(x/r2);
    while (1){
        if(fabs(angle2+2*pi-angle)<fabs(angle2-angle)){
            angle2+=2*pi;
            continue;
        }else if(fabs(angle2-2*pi-angle)<fabs(angle2-angle)){
            angle2-=2*pi;
            continue;
        }
        break;
    } 
    angle=angle2;
    dots[i]={asin(z/r),angle2};
    }
}
void process_value_lr(vector<vector<double>> &dots, vector<double> centre,double r){
    static double angle=0;
    for (size_t i = 0; i < dots.size(); i++)
    {
        if(dots[i].size()<3){
            cout<<"invalid input"<<endl;
            throw "invalid input";
        }
    double x=dots[i][0]-centre[0],y=dots[i][1]-centre[1],z=dots[i][2]-centre[2];
    double r2=sqrt(x*x+z*z);
    if(r2==0){
        dots[i]= {asin(y/r),angle};
        continue;
    }
    double angle2=z>0 ? acos(x/r2):2*pi-acos(x/r2);
    while (1){
        if(fabs(angle2+2*pi-angle)<fabs(angle2-angle)){
            angle2+=2*pi;
            continue;
        }else if(fabs(angle2-2*pi-angle)<fabs(angle2-angle)){
            angle2-=2*pi;
            continue;
        }
        break;
    } 
    angle=angle2;
    dots[i]={asin(y/r),angle2};
    }
        
}

vector<double> process_derivation_ud(vector<double> dots,vector<double> values, int order,double r,double rate){
    double coeff_11,coeff_21,coeff_22,coeff_31,coeff_32;
    int sign=order%4<2 ? 1:-1;
    if(order%2==1) {
        coeff_11=cos(dots[0])*r*sign;
        coeff_32=cos(dots[1])*r*sign*cos(dots[0]);
    }
    else {
        coeff_11=sin(dots[0])*r*sign;
        coeff_32=sin(dots[1])*r*sign*cos(dots[0]);
    }
    if(order%4==1||order%4==2) sign=-1;
    else sign=1;
    if(order%2==0){
        coeff_21=cos(dots[0])*sign*r*cos(dots[1]);
        coeff_31=cos(dots[0])*sign*r*sin(dots[1]);
        coeff_22=cos(dots[1])*r*sign*cos(dots[0]);
    }else{
        coeff_21=sin(dots[0])*sign*r*cos(dots[1]);
        coeff_31=sin(dots[0])*sign*r*sin(dots[1]);
        coeff_22=sin(dots[1])*r*sign*cos(dots[0]);
    }
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(3, 2);
    Eigen::VectorXd target(3);
    matrix(0,0)=coeff_11;
    matrix(1,0)=coeff_21;
    matrix(1,1)=coeff_22;
    matrix(2,0)=coeff_31;
    matrix(2,1)=coeff_32;
    target[0]=values[0]*rate;
    target[1]=values[1]*rate;
    target[2]=values[2]*rate;
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
    if(lu_decomp.rank()<2){
        cout<<"incorrect derivative input"<<endl;
        throw "incorrect derivative input";
    }
    Eigen::VectorXd solution = matrix.colPivHouseholderQr().solve(target);
    vector<double> result(solution.data(), solution.data() + solution.size());
    return result;
}

vector<double> process_derivation_lr(vector<double> dots,vector<double> values, int order,double r,double rate){
    double coeff_11,coeff_21,coeff_22,coeff_31,coeff_32;
    int sign=order%4<2 ? 1:-1;
    if(order%2==1) {
        coeff_11=cos(dots[0])*r*sign;
        coeff_32=cos(dots[1])*r*sign*cos(dots[0]);
    }
    else {
        coeff_11=sin(dots[0])*r*sign;
        coeff_32=sin(dots[1])*r*sign*cos(dots[0]);
    }
    if(order%4==1||order%4==2) sign=-1;
    else sign=1;
    if(order%2==0){
        coeff_21=cos(dots[0])*sign*r*cos(dots[1]);
        coeff_31=cos(dots[0])*sign*r*sin(dots[1]);
        coeff_22=cos(dots[1])*r*sign*cos(dots[0]);
    }else{
        coeff_21=sin(dots[0])*sign*r*cos(dots[1]);
        coeff_31=sin(dots[0])*sign*r*sin(dots[1]);
        coeff_22=sin(dots[1])*r*sign*cos(dots[0]);
    }
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(3, 2);
    Eigen::VectorXd target(3);
    matrix(0,0)=coeff_11;
    matrix(1,0)=coeff_21;
    matrix(1,1)=coeff_22;
    matrix(2,0)=coeff_31;
    matrix(2,1)=coeff_32;
    target[0]=values[0]*rate;
    target[1]=values[2]*rate;
    target[2]=values[1]*rate;
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
    if(lu_decomp.rank()<2){
        cout<<"incorrect derivative input"<<endl;
        throw "incorrect derivative input";
    }
    Eigen::VectorXd solution = matrix.colPivHouseholderQr().solve(target);
    vector<double> result(solution.data(), solution.data() + solution.size());
    return result;
}