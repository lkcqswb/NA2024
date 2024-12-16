#include "../problems_1_7/problem3/BSpline_2_3_complete.hpp"
#include "../include/json.hpp"
#include <vector>
using namespace std;
using json = nlohmann::json;
double f(double x){
    return (double)(1)/(1+x*x);
}


int main(){
    vector<double> knots1={-5,-4,-3,-2,-1,0,1,2,3,4,5},knots2={-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5};
    vector<double> value1,value2;
    for (size_t i = 0; i < knots1.size(); i++) value1.push_back(f(knots1[i]));
    for (size_t i = 0; i < knots2.size(); i++) value2.push_back(f(knots2[i]));  
    json j1={
            {"data points", knots1},
            {"values",value1},
            {"diff_values",{double(10)/(26*26),-double(10)/(26*26)}},
            {"range", {
                {"end", 5},
                {"begin", -5}
            }}
        };
    json j2={
            {"data points", knots2},
            {"values",value2},
            {"boundary_values",{f(-5),f(5)}},
            {"range", {
                {"end", 5},
                {"begin", -5}
            }}
        };
    
    BSpline_uniform_23 bs23(j1);
    BSpline_uniform_12 bs12(j2);

    ofstream outfile("Pc.txt", ios::trunc);
    for (double t = -5; t <= 5; t+=0.1) outfile << t << " " << bs23.get_value(t) << endl;
    outfile << "#END# " <<"B^2_3"<< endl;
    for (double t = -4.5; t <= 4.5; t+=0.1) outfile << t << " " << bs12.get_value(t) << endl;
    outfile << "#END# " <<"B^1_2"<< endl;
    for (double t = -5; t <= 5; t+=0.1) outfile << t << " " << f(t) << endl;
    outfile << "#END#The function itself"<< endl;
    string command = "python plot.py Pc.txt";
    system(command.c_str());
    
}