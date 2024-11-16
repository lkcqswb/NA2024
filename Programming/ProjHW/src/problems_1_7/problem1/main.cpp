#include "ppform_s_1_0.hpp"
#include "BSpline_s_1_0.hpp"
#include <iostream>
#include <fstream>
#include "../../include/json.hpp"
using json = nlohmann::json;


int main() {
    json j;
    std::ifstream file("config.json"); 
    if (!file.is_open()) {
            std::cerr << "无法打开文件" << std::endl;
            return 1;
    }
    j = json::parse(file);
    file.close();
    ppform_1_0 pp(j);
    BSpline_1_0 bs(j);

    double t;
    ofstream outfile1("ppform_s_1_0.txt", ios::trunc);
    ofstream outfile2("BSpline_s_1_0.txt", ios::trunc);

    for (t = j["range"]["begin"]; t <= j["range"]["end"]; t+=0.01) {
            outfile1 << t << " " << pp.get_value(t) << endl;
            outfile2 << t << " " << bs.get_value(t) << endl;
    }



    string command = "python plot.py ppform_s_1_0.txt";
    system(command.c_str());
    command = "python plot.py BSpline_s_1_0.txt";
    system(command.c_str());


    return 0;

}