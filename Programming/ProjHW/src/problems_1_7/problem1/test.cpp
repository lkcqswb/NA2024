#include "../../ppForm_a_BSpline/BSpline.hpp"
#include "../../ppForm_a_BSpline/ppForm.hpp"
#include "../../include/json.hpp"
#include <fstream>
#include <iostream>

using json = nlohmann::json;
int main(){
    std::ifstream file("config.json");
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << "config.json" << std::endl;
        return 0;
    }

    json j;
    file >> j;

    ppForm bs(j);
    double eps=1e-6;
    
    vector<double> knots=j["data points"][0].get<vector<double>>();
    for (size_t i = 1; i < knots.size()-1; i++)
    {
        cout<<(-bs.get_value(knots[i]-eps)[0]+bs.get_value(knots[i])[0])/eps<<endl;
        cout<<(bs.get_value(knots[i]+eps)[0]-bs.get_value(knots[i])[0])/eps<<endl;
    }
    for (size_t i = 1; i < knots.size()-1; i++)
    {
        cout<<(bs.get_value(knots[i]-2*eps)[0]-2*bs.get_value(knots[i]-eps)[0]+bs.get_value(knots[i])[0])/(eps*eps)<<endl;
        cout<<(bs.get_value(knots[i]+2*eps)[0]-2*bs.get_value(knots[i]+eps)[0]+bs.get_value(knots[i])[0])/(eps*eps)<<endl;
    }
    cout<<(bs.get_value(knots[3]-2*eps)[0]-2*bs.get_value(knots[3]-eps)[0]+bs.get_value(knots[3])[0])/(eps*eps)<<endl;
    cout<<(bs.get_value(knots[2]-2*eps)[0]-2*bs.get_value(knots[2]-eps)[0]+bs.get_value(knots[2])[0])/(eps*eps)<<endl;
    cout<<(bs.get_value(knots[1]-2*eps)[0]-2*bs.get_value(knots[1]-eps)[0]+bs.get_value(knots[1])[0])/(eps*eps)<<endl;
    cout<<(bs.get_value(knots[0]+3*eps)[0]-2*bs.get_value(knots[0]+2*eps)[0]+bs.get_value(knots[0]+eps)[0])/(eps*eps)<<endl;
    
   //cin>>t;
   //cout<<bs.get_value(t)[0];
   
}