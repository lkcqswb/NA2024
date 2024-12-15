#include "..\problems_1_7\problem2\ppform_2_3_natural.hpp"
#include "..\problems_1_7\problem3\BSpline_2_3_natural.hpp"
#include <iostream>
#include <vector>
#include <fstream> 
#include <cmath>  
#include <cstdlib> 
using namespace std;

vector<double> construct(int N, double start, double end) {
    vector<double> knots;
    for (int i = 0; i < N; i++) {
        knots.push_back((double)i / (N - 1) * (end - start) + start);
    }
    return knots;
}

double f(double x) {
    return 1 / (1 + 25 * x * x);//构建函数值表格
}

vector<double> construct_f(vector<double> knots) {
    vector<double> values;
    for (size_t i = 0; i < knots.size(); i++) {
        values.push_back(f(knots[i]));
    }
    return values;
}

double test_max_norm_p(int N, double start, double end) {
    vector<double> knots = construct(N, start, end);
    vector<double> f_values = construct_f(knots);
    json j = {
        {"dimension", 1},
        {"order", 3},
        {"boundary condition", "natural"},
        {"data points", knots},
        {"function values", f_values},
        {"range", {
            {"end", 1},
            {"begin", -1}
        }}
    };
    ppform_2_3_natural pp(j);
    double max_error = 0;
    for (size_t i = 0; i < knots.size() - 1; i++) {
        double result = pp.get_value((knots[i] + knots[i + 1]) / 2);
        double result2 = f((knots[i] + knots[i + 1]) / 2);
        if (abs(result2 - result) > max_error) max_error = abs(result2 - result);
    }
    return max_error;
}

double test_max_norm_B(int N, double start, double end) {
    vector<double> knots = construct(N, start, end);
    vector<double> f_values = construct_f(knots);
    json j = {
        {"dimension", 1},
        {"order", 3},
        {"boundary condition", "natural"},
        {"data points", knots},
        {"function values", f_values},
        {"range", {
            {"end", 1},
            {"begin", -1}
        }}
    };
    BSpline_2_3_natural bs(j);
    double max_error = 0;
    for (size_t i = 0; i < knots.size() - 1; i++) {
        double result = bs.get_value((knots[i] + knots[i + 1]) / 2);
        double result2 = f((knots[i] + knots[i + 1]) / 2);
        if (abs(result2 - result) > max_error) max_error = abs(result2 - result);
    }
    return max_error;
}

void write_errors_to_file() {
    int n[5] = {6, 11, 21, 41, 81};
    ofstream outfile("errors.csv"); 
    outfile << "N,Error_PPForm,Error_BSpline\n"; 

    // 计算误差并写入文件
    for (size_t i = 0; i < 5; i++) {
        double error_pp = test_max_norm_p(n[i], -1, 1);
        double error_bs = test_max_norm_B(n[i], -1, 1);
        outfile << n[i] << "," << error_pp << "," << error_bs << endl;
    }

    outfile.close();
    cout << "Error data has been written to errors.csv\n";
}

void call_python_script() {
    const char* pythonScript = "plot_errors.py";
    system(("python " + string(pythonScript)).c_str()); 
}

int main() {
    write_errors_to_file();  
    call_python_script();   
    return 0;
}