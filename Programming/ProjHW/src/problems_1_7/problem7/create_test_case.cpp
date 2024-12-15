#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../../include/json.hpp"
using namespace std;
using json = nlohmann::json;
#define pi 3.1415926535
std::vector<std::vector<double>> generateSphericalSpiral(int num_points, double height, double radius) {
    std::vector<std::vector<double>> points;
    
    for (size_t i = 0; i < num_points; ++i) {
        double theta = i * (2 * pi / 15); 
        double phi = (pi / height) * i; 

        double x = radius * sin(phi) * cos(theta);
        double y = radius * sin(phi) * sin(theta);
        double z = radius * cos(phi);
        
        points.push_back({x, y, z});
    }
    
    return points;
}


int main() {
    vector<vector<double>> points;
    double radius=1;
    vector<double> centre={0,0,0};
    //////
    //for (size_t i = 0; i < 20; i++)
    //{
    //    points.push_back({cos(i),sin(i)*cos(2*i-1),sin(i)*sin(2*i-1)});
    //}
    points = generateSphericalSpiral(100, 30, radius);

    json j;
    j["order"] = 3;
    j["radius"]= radius,
    j["centre"] = centre;
    j["boundary condition"]={
        {"exists",{1,points.size()-2}}
    };
    j["points"]=points;
    j["range"]={
        {"begin",0},
        {"end",1}
    };
    std::ofstream file("config.json");
    if (file.is_open()) {
        file << j.dump(4); 
        file.close();
        std::cout << "JSON data written to config.json\n";
    } else {
        std::cerr << "Unable to open file for writing\n";
    }

    return 0;
}