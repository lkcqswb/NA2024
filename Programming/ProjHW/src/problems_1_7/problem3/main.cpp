#include "BSpline_2_3_periodic.hpp"
#include "BSpline_2_3_complete.hpp"
#include "BSpline_2_3_natural.hpp"
#include "BSpline_2_3_not_a_knot.hpp"
#include "BSpline_2_3_specific_2_deri.hpp"
#include<iostream>
#include <vector>
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
    double start,end;
    vector<double> knots;
    if (!j["data points"].is_null()){
        knots=j["data points"].get<vector<double>>();
        if(knots.empty()||knots.size()<2) {cout<<"enadequate knots"<<endl; throw "enadequate datapoints";}
    }else {cout<<"no datapoints"<<endl; throw "no datapoints";}

    if (!j["range"]["begin"].is_null()) start=j["range"]["begin"];
    else start=knots[0];
    if (!j["range"]["end"].is_null()) end=j["range"]["end"];
    else end=knots[knots.size()-1];
    


    ofstream outfile1("BSpline_2_3.txt", ios::trunc);
    double t;
    if(!j["boundary condition"].is_null()){
        if(j["boundary condition"]=="complete"){
            BSpline_2_3_complete pp(j);
            for (t = start; t <= end; t+=0.01) {
                outfile1 << t << " " << pp.get_value(t) << endl;
            }
        
        }else if(j["boundary condition"]=="natural"){
                BSpline_2_3_natural pp(j);
                for (t = start; t <= end; t+=0.01) {
                    outfile1 << t << " " << pp.get_value(t) << endl;
                }
            
        }else if(j["boundary condition"]=="periodic"){
                BSpline_2_3_periodic pp(j);
                for (t = start; t <= end; t+=0.01) {
                    outfile1 << t << " " << pp.get_value(t) << endl;
                }
            
        }else if(j["boundary condition"]=="not_a_knot"){
                BSpline_2_3_not_a_knot pp(j);
                for (t = start; t <= end; t+=0.01) {
                    outfile1 << t << " " << pp.get_value(t) << endl;
                }
            
        }else if(j["boundary condition"]=="specific_2_deri"){
                BSpline_2_3_specific_2_deri pp(j);
                for (t = start; t <= end; t+=0.01) {
                    outfile1 << t << " " << pp.get_value(t) << endl;
                }
            
        }
    }



    string command = "python plot.py BSpline_2_3.txt";
    system(command.c_str());



    return 0;

}