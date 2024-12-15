#include "BSpline_s_1_0.hpp"
#include "ppform_s_1_0.hpp"
#include<vector>
#include<fstream>

using json = nlohmann::json;
using namespace std;
int main() {
    json j;
    std::ifstream file("config.json"); 
    if (!file.is_open()) {
            std::cerr << "无法打开文件" << std::endl;
            return 1;
    }
    j = json::parse(file);
    file.close();
    ppForm pp(j,2);
    BSpline bs(j);

    double t;
    ofstream outfile1("p_s_1_0.txt", ios::trunc);
    ofstream outfile2("B_s_1_0.txt", ios::trunc);

    for (t = j["range"]["begin"]; t <= j["range"]["end"]; t+=0.1) {
        vector<double> result1=pp.get_value(t);
        vector<double> result2=bs.get_value(t);
        outfile1 <<t<<" "<<result1[0]  << endl;
        outfile2 <<t<<" "<<result2[0]  << endl;
    }



    string command = "python plot.py p_s_1_0.txt";
    system(command.c_str());
    command = "python plot.py B_s_1_0.txt";
    system(command.c_str());
}