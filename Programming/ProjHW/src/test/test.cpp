#include "../include/json.hpp"
#include "../ppForm_a_BSpline/ppForm.hpp"
#include "../ppForm_a_BSpline/BSpline.hpp"
#include<vector>
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
    ofstream outfile1("sphere_fitting_p.txt", ios::trunc);
    ofstream outfile2("sphere_fitting_B.txt", ios::trunc);

    for (t = j["range"]["begin"]; t <= j["range"]["end"]; t+=0.1) {
        vector<double> result1=pp.get_value(t);
        vector<double> result2=bs.get_value(t);
        cout << result1[0] << endl;
        cout << result2[0] << endl;
    }



    string command = "python plot.py sphere_fitting_p.txt";
    system(command.c_str());
    command = "python plot.py sphere_fitting_B.txt";
    system(command.c_str());
}