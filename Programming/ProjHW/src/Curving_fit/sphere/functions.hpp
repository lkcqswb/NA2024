#pragma once
#include<vector>
#include<iostream>
#include"../../include/Eigen/Dense"
#define pi 3.1415926535
using namespace std;
using namespace Eigen;
#include <cstdlib>
#include <ctime>
#include <random>
#include <cmath>

// 计算旋转矩阵，使得 point 旋转到北极
Matrix3d get_rotation_matrix_to_north_pole(const Vector3d &point) {
    Vector3d north_pole(0, 0, 1);

    Vector3d rotation_axis = point.cross(north_pole);

    if (rotation_axis.norm() < 1e-6) {
        return Matrix3d::Identity(); // 返回单位矩阵
    }

    rotation_axis.normalize();


    double angle = acos(point.dot(north_pole) / (point.norm() * north_pole.norm()));

 
    AngleAxisd angle_axis(angle, rotation_axis);
    return angle_axis.toRotationMatrix(); // 返回正交矩阵
}

// 检查点是否在 dots 中
bool is_point_in_dots(const Vector3d &point, const vector<vector<double>> &dots) {
    for (const auto &dot : dots) {
        if (dot.size() < 3) continue;
        Vector3d d(dot[0], dot[1], dot[2]);
        if ((d - point).norm() < 1e-6) return true;
    }
    return false;
}

// 随机生成单位球面上的点
Vector3d generate_random_point_on_sphere() {
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_real_distribution<double> dis(0.0, 1.0);

    double u = dis(gen);
    double v = dis(gen);
    double theta = 2 * pi * u; // 均匀分布的角度
    double phi = acos(2 * v - 1); // 均匀分布的极角

    double x = sin(phi) * cos(theta);
    double y = sin(phi) * sin(theta);
    double z = cos(phi);
    return Vector3d(x, y, z);
}

Matrix3d process_value(vector<vector<double>> &dots, vector<double> centre, double r) {
    Vector3d random_point, symmetric_point,generate_rotate;
    // 不断随机选择点，直到满足条件
    while (true) {
        random_point = generate_random_point_on_sphere();
        generate_rotate=random_point;
        symmetric_point = -random_point;

        random_point[0]=centre[0]+random_point[0]*r;
        random_point[1]=centre[1]+random_point[1]*r;
        random_point[2]=centre[2]+random_point[2]*r;
        symmetric_point[0]=centre[0]+symmetric_point[0]*r;
        symmetric_point[1]=centre[1]+symmetric_point[1]*r;
        symmetric_point[2]=centre[2]+symmetric_point[2]*r;
        // 检查随机点及其对称点是否在 dots 中
        if (!is_point_in_dots(random_point, dots) && !is_point_in_dots(symmetric_point, dots)) {
            break; // 找到了合适的点，退出循环
        }
    }
    Matrix3d rotation_matrix = get_rotation_matrix_to_north_pole(generate_rotate);
    for (size_t i = 0; i < dots.size(); i++) {
        Vector3d point(dots[i][0] - centre[0], dots[i][1] - centre[1], dots[i][2] - centre[2]);
        Vector3d rotated_point = rotation_matrix * point;
        dots[i] = {rotated_point[0] + centre[0], rotated_point[1] + centre[1], rotated_point[2] + centre[2]};
    }
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
    return rotation_matrix;
}


double sign(double x){
    if(x>=0) return 1;
    return -1;
}


vector<double> process_derivation(vector<double> dots,vector<double> values, int order,double r,double rate,Matrix3d rotation_matrix){
    double coeff_11,coeff_21,coeff_22,coeff_31,coeff_32;
    Vector3d point(values[0], values[1] , values[2]);
    Vector3d rotated_point = rotation_matrix * point;
    values={rotated_point[0],rotated_point[1],rotated_point[2]};

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
