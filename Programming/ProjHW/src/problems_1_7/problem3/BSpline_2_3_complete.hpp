#pragma once

#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/BSpline.hpp"
#include<vector>
using json = nlohmann::json;
using namespace std;




class BSpline_2_3_complete: public BSpline{
private:
    vector<double> knots; 
    double start,end;
    vector<int> dir,dot;
    
    vector<int> generate_seq(int n){
        vector<int> result;
        for (size_t i = 0; i < n; i++)
        {
            result.push_back(i);
        }
        result.push_back(0);
        result.push_back(n-1);//对应两个2阶导的条件中点序号
        return result;
    }
    vector<int> generate_orders(int n){
        vector<int> result;
        for (size_t i = 0; i < n; i++)
        {
            result.push_back(0);
        }
        result.push_back(1);
        result.push_back(1);//对应两个2阶导的条件中导数阶
        return result;
    }

    vector<vector<double>> generate_func_value(vector<double> f_values,vector<double> derivations){
        if(derivations.size()<2){
            cout<<"enadequate derivation information"<<endl;
            throw "enadequate derivation information";
        }
        vector<vector<double>> result;
        for (size_t i = 0; i < f_values.size(); i++)
        {
            result.push_back({f_values[i]});
        }
        result.push_back({derivations[0]});
        result.push_back({derivations[1]});//表示最后两个一阶导值。
        return result;
    }
public:
    BSpline_2_3_complete(vector<double>input_knots,vector<double> f_values,vector<double> derivations,double istart,double iend):
    BSpline({
            {"dimension", 1},
            {"order", 3},
            {"boundary condition", {
                {"values", {
                    generate_seq(input_knots.size()),
                    generate_orders(input_knots.size()),
                    generate_func_value(f_values,derivations)
                }}
            }},
            {"data points", input_knots},
            {"range", {
                {"end", iend},
                {"begin", istart}
            }}
        }
    ),knots(input_knots),start(istart),end(iend){};
    double get_value(double t){
        return BSpline::get_value(t)[0];
    }
};




