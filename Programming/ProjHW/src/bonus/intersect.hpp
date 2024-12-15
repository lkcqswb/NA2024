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

// 计算点积
double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

// 计算向量减法
std::vector<double> subtract(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// 判断两个向量是否相等，允许误差
bool equal(const std::vector<double>& a, const std::vector<double>& b, double tolerance = 1e-6) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i] - b[i]) > tolerance) return false;
    }
    return true;
}

// 判断两条线段是否相交（适用于 N 维空间）
bool segments_intersect_nd(const PointND& p1, const PointND& q1, const PointND& p2, const PointND& q2, double epsilon = 1e-9) {
    auto d1 = subtract(q1.coords, p1.coords);
    auto d2 = subtract(q2.coords, p2.coords);
    auto r = subtract(p1.coords, p2.coords);

    double a = dot(d1, d1);
    double b = dot(d1, d2);
    double c = dot(d2, d2);
    double d = dot(d1, r);
    double e = dot(d2, r);

    double denom = a * c - b * b;

    if (std::fabs(denom) > epsilon) { // 不平行的情况
        double t = (b * e - c * d) / denom;
        double s = (a * e - b * d) / denom;

        // 检查参数 t 和 s 是否在 [0, 1] 范围内
        if (t >= 0 && t <= 1 && s >= 0 && s <= 1) return true;
    } else { // 平行的情况
        auto r_proj = dot(d1, r);
        auto q_proj = dot(d1, subtract(q2.coords, p1.coords));
        if (std::fabs(r_proj) < epsilon && std::fabs(q_proj) < epsilon) { // 共线，检查投影区间是否重叠
            double d1_len = dot(d1, d1);
            double d2_len = dot(d2, d2);
            return (d >= 0 && d <= d1_len) || (e >= 0 && e <= d2_len);
        }
    }

    return false;
}
template <typename T>
T clamp(T value, T low, T high) {
    return (value < low) ? low : (value > high) ? high : value;
}
// 判定是否存在自相交
bool detect_intersect(const json& j) {
    // 初始化样条曲线
    curve_fitting_B B(j, "Chord");
    double start = j["range"].contains("begin") ? j["range"]["begin"].get<double>() : 0;
    double end = j["range"].contains("end") ? j["range"]["end"].get<double>() : 1;

    int dimension = j["dimension"].get<int>();
    int num_samples = 2000; // 增加样本点数量，提高检测精度

    // 离散化样条曲线
    std::vector<PointND> points;
    for (int i = 0; i <= num_samples; ++i) {
        double t = start + i * (end - start) / num_samples;
        t = clamp(t, start, end);
        auto value = B.get_value(t);
        points.push_back({value});
    }

    // 检查首尾点重合
    if (equal(points.front().coords, points.back().coords)) {
        return true; // 存在自相交
    }

    // 检查线段是否相交
    int n = points.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 2; j < n - 1; ++j) { // 避免检查相邻线段
            if (segments_intersect_nd(points[i], points[i + 1], points[j], points[j + 1])) {
                return true; // 存在自相交
            }
        }
    }

    // 检查首尾线段是否相交
    if (segments_intersect_nd(points[0], points[1], points[n - 2], points[n - 1])) {
        return true;
    }

    return false; // 没有自相交
}


