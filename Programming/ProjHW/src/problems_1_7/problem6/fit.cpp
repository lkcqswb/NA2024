#include "../../Curving_fit/plane/curve_fitting_B.hpp"
#include "../../Curving_fit/plane/curve_fitting_p.hpp"
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
    std::string fileContent((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    try {
        j = json::parse(fileContent);
    } catch (const json::parse_error& e) {
        std::cerr << "JSON 解析错误: " << e.what() << std::endl;
        return 1; 
    }
    file.close();
    
    curve_fitting_B bs(j);
    curve_fitting_p pp(j);
    
    double t;
    ofstream outfile1("curve_fitting_p.txt", ios::trunc);
    ofstream outfile2("curve_fitting_B.txt", ios::trunc);

    for (t = j["range"]["begin"]; t <= j["range"]["end"]; t+=0.01) {
        vector<double> result1=pp.get_value(t);
        vector<double> result2=bs.get_value(t);
        outfile1 << result1[0] << " " << result1[1] << endl;
        outfile2 << result2[0] << " " << result2[1] << endl;
    }



    string command = "python plot.py curve_fitting_p.txt";
    system(command.c_str());
    command = "python plot.py curve_fitting_B.txt";
    system(command.c_str());


    return 0;

}