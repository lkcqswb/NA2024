#pragma once
#include<vector>
#include<iostream>
#include"../../include/Eigen/Dense"
#define pi 3.1415926535
using namespace std;


vector<double> process_value(vector<double> dots, vector<double> centre,double r){
    static double angle=0;
    if(dots.size()<3){
        cout<<"invalid input"<<endl;
        throw "invalid input";
    }
    double x=dots[0]-centre[0],y=dots[1]-centre[1],z=dots[2]-centre[2];
    double r2=sqrt(x*x+y*y);
    if(r2==0) return {asin(z/r),angle};
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
    return {asin(z/r),angle2};
}

vector<double> process_derivation(vector<double> dots,vector<double> values, int order,double r){
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
    target=Eigen::Map<const Eigen::VectorXd>(values.data(),3);
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
    if(lu_decomp.rank()<2){
        cout<<"incorrect derivative input"<<endl;
        throw "incorrect derivative input";
    }
    Eigen::VectorXd solution = matrix.colPivHouseholderQr().solve(target);
    vector<double> result(solution.data(), solution.data() + solution.size());
    return result;
}