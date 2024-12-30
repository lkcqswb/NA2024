#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include "../ppForm_a_BSpline/BSpline.hpp"  
#include "../ppForm_a_BSpline/ppForm.hpp"
#include "../include/json.hpp"
using namespace std;



using Function = std::function<double(double)>;

// 计算最大误差
double compute_max_error(const vector<double>& true_values, const vector<double>& approx_values) {
    double max_error = 0.0;
    for (size_t i = 0; i < true_values.size(); ++i) {
        max_error = max(max_error, abs(true_values[i] - approx_values[i]));
    }
    return max_error;
}


void convergence_analysis(const Function& true_function, const pair<double, double>& domain, const vector<int>& n_values,string kinds ) {
    cout << setw(10) << "n" 
         << setw(20) << "Max Error"
         << setw(20) << "Convergence Order" << endl;
    cout << string(50, '-') << endl;

    double prev_error = -1; // 用于存储上一个步长的误差
    double prev_h = -1;     // 用于存储上一个步长的 h 值

    // 遍历节点
    for (int n : n_values) {
        vector<double> x(n),orders(n),num(n);
        vector<vector<double>> y(n);
        // 均匀生成节点
        double h = (domain.second - domain.first) / (n - 1);
        for (int i = 0; i < n; ++i) {
            x[i] = domain.first + i * h;
            orders[i]=0;
            num[i]=i;
            y[i] = {true_function(x[i])};
        }

        json j={
            {"dimension",1},
            {"order", 3},
            {"boundary condition", {
                {"values", {
                    num,
                    orders,
                    y
                }},
                {"exists",{1,n-2}}
            }},
            {"data points", x},
            {"range", {
                {"end", domain.second},
                {"begin", domain.first}
            }}
        };
        BSpline B(j);
        ppForm p(j);

        // 计算误差
        vector<double> fine_x, true_y, approx_y;
        int fine_n = 1000;
        double fine_h = (domain.second - domain.first) / (fine_n - 1);
        for (int i = 0; i < fine_n; ++i) {
            double x_val = domain.first + i * fine_h;
            fine_x.push_back(x_val);
            true_y.push_back(true_function(x_val)); 
            if(kinds=="B") approx_y.push_back(B.get_value(x_val)[0]);    
            if(kinds=="P") approx_y.push_back(p.get_value(x_val)[0]);  
        }

        // 计算最大误差
        double max_error = compute_max_error(true_y, approx_y);

        // 计算收敛阶
        double order = -1; 
        if (prev_error > 0 && prev_h > 0) {
            order = log(prev_error / max_error) / log(prev_h / h);
        }

        // 输出步长、误差和收敛阶
        cout << setw(10) << n
             << setw(20) << max_error
             << setw(20) << (order > 0 ? to_string(order) : "-") << endl;

        prev_error = max_error;
        prev_h = h;
    }
}


int main() {
    //y = sin(x)
    Function true_function = [](double x) { return sin(x); };

    // 定义域
    pair<double, double> domain = {0.0, M_PI};
    vector<int> n_values;
    // 节点数
    for (size_t i = 2; i < 8; i++)
    {
        n_values.push_back(pow(2,i));
    }

    convergence_analysis(true_function, domain, n_values,"B");
    convergence_analysis(true_function, domain, n_values,"P");

    return 0;
}