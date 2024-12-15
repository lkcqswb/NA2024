#include "../Curving_fit/plane/curve_fitting_B.hpp"
#include "../Curving_fit/plane/curve_fitting_p.hpp"
#pragma once
/*
json 格式为{
    "dimension":2,
    "order": 3,
    "boundary condition": {
        "equals": [],
        "values": [],
        "exists": [1,2]
    },
    "points": 
        [[0,2],[ 2,7],[ 4,4],[9,10]],
    "range":{
        "begin":0,
        "end":1
    }
}(仅为样例，维数等信息可随意更改)
*/
#include <vector>
#include <cmath>
#include <limits>

struct PointND {
    std::vector<double> coords;
};


double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

// 向量减法
std::vector<double> subtract(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

bool equal(vector<double> a,vector<double> b){
    for (size_t i = 0; i < a.size(); i++)
    {
        if(fabs(a[i]-b[i])>1e-6) return false;
    }
    return true;
}

// 判断两条线段是否相交
bool segments_intersect_nd(const PointND& p1, const PointND& q1, const PointND& p2, const PointND& q2) {
    const double epsilon = 1e-9;


    auto d1 = subtract(q1.coords, p1.coords);
    auto d2 = subtract(q2.coords, p2.coords);
    auto r = subtract(p1.coords, p2.coords);

  
    double a = dot(d1, d1);
    double b = dot(d1, d2);
    double c = dot(d2, d2);
    double d = dot(d1, r);
    double e = dot(d2, r);


    double denom = a * c - b * b;
    double t, s;

    if (fabs(denom) > epsilon) { // 不平行
        t = (b * e - c * d) / denom;
        s = (a * e - b * d) / denom;

       
        if (t < 0 || t > 1 || s < 0 || s > 1) return false;
    } else { // 平行
      
        double d1_proj = dot(d1, r);
        double d2_proj = dot(d2, subtract(q1.coords, p2.coords));
        if (fabs(d1_proj) > epsilon || fabs(d2_proj) > epsilon) return false;
    }

    // 距离小于阈值，认为相交
    auto p1_proj = subtract(p1.coords, subtract(d1, d2));
    double distance = dot(p1_proj, p1_proj);
    return distance < epsilon;
}


bool detect_intersect(json j) {
    curve_fitting_B B(j, "Chord");
    double start,end;
    if (!j["range"]["begin"].is_null()) start=j["range"]["begin"];
        else start=0;
    if (!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=1;
    int dimension = j["dimension"];

    int num_samples = 1000; // 增加离散点数
    std::vector<PointND> points;

    // 离散化样条曲线
    for (int i = 0; i <= num_samples; ++i) {
        double t = std::max(std::min(start + i * (end - start) / num_samples, end), start);
        auto value = B.get_value(t);
        points.push_back({value});
    }

    // 检查点重合（头尾相连情况）
    if (equal(points.front().coords, points.back().coords)) {
        return true; // 存在自相交
    }

    // 检查线段交叉
    int n = points.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 2; j < n - 1; ++j) { // 避免相邻线段检查
            if (segments_intersect_nd(points[i], points[i + 1], points[j], points[j + 1])) {
                return true; // 存在自相交
            }
        }
    }

    // 检查头尾线段交叉
    if (segments_intersect_nd(points[0], points[1], points[n - 2], points[n - 1])) {
        return true; // 存在自相交
    }

    return false; // 没有自相交
}

