#ifndef PPFORM
#define PPFORM
#include<iostream>
#include "../include/json.hpp"
#include "ppForm_options/solve_boundary_condition.hpp"
#include<fstream>
#include<string>
#include<vector>
#include <algorithm>
using namespace std;
using json = nlohmann::json;

class ppForm
{
private:
    vector<double> knots;
    vector<vector<vector<double>>> coef;
    double begin,end;
    int order;
    int dimension;
public:
    ppForm(json j);
    ~ppForm();
    vector<double> get_value(double x);
};

ppForm::ppForm(json j)
{
    //检查输入的json是否合法
    if (!j["dimension"].is_null()) dimension=j["dimension"];
    if (!j["order"].is_null()) order=j["order"];
    vector<vector<double>> f_values;
    if (!j["data points"].is_null()){
        knots=j["data points"][0].get<vector<double>>();
        if(knots.empty()||knots.size()<2) {cout<<"enadequate knots"<<endl; throw "enadequate datapoints";}
        f_values=j["data points"][1].get<vector<vector<double>>>();
        if(knots.size()!=f_values.size()) {cout<<"illegal shape"<<endl; throw "illegal shape";}
        for(auto &item:f_values){if((int)item.size()<dimension){cout<<"illegal shape"<<endl; throw "illegal shape";}}//如果函数值不符合维数要求
    }else {cout<<"no datapoints"<<endl; throw "no datapoints";}
    //对结点排序
    vector<pair<double, vector<double>>> pairs;
    for (size_t i = 0; i < knots.size(); ++i) {
        pairs.push_back(make_pair(knots[i], f_values[i]));
    }
    sort(pairs.begin(), pairs.end(), [](const std::pair<double, vector<double>>& a, const std::pair<double, vector<double>>& b) {
        return a.first < b.first;
        });
    for (size_t i = 0; i < pairs.size(); ++i) {
        knots[i] = pairs[i].first;
        f_values[i] = pairs[i].second;
    }
    //排序结束



    if (!j["range"]["begin"].is_null()) begin=j["range"]["begin"];
    else begin=knots[0];
    if (!j["range"]["end"].is_null()) end=j["range"]["end"];
    else end=knots[knots.size()-1];
    //结束检查
    
    for(int cur_dimen=0;cur_dimen<dimension;cur_dimen++){
        vector<double> f_value_dimension;
        for(size_t j=0;j<knots.size();j++) f_value_dimension.push_back(f_values[j][cur_dimen]);//该维度的所有函数值
             
        vector<int> dots1={},dots2={},difforder1={},difforder2={};
        
        vector<int> dots={},difforder={};
        vector<double> values_dimension={};

        vector <int> exist_dot={};
        

        if(!j["boundary condition"]["equals"].is_null()){
            json eq = j["boundary condition"]["equals"];
            if (!eq.empty()){
                if(eq.size()<4){
                    cout<<"illegal boundary condition"<<endl;
                    throw "illegal boundary condition";
                }
                dots1=eq[0].get<vector<int>>();
                difforder1=eq[1].get<vector<int>>();
                dots2=eq[2].get<vector<int>>();
                difforder2=eq[3].get<vector<int>>();
                if(dots1.size()!=difforder1.size()||dots2.size()!=difforder2.size()||dots1.size()!=dots2.size()){
                    cout<<"The number of left endpoint sequences, derivative orders, right endpoint sequences, and derivative orders should be the same."<<endl;
                    throw "invalid shapes";
                }
            }

            
        }
        if(!j["boundary condition"]["values"].is_null()){
            json va = j["boundary condition"]["values"];
            vector<vector<double>> value={};
            if (!va.empty()){
                if(va.size()<3){
                    cout<<"illegal boundary condition"<<endl;
                    throw "illegal boundary condition";
                }
                dots=va[0].get<vector<int>>();
                difforder=va[1].get<vector<int>>();
                value=va[2].get<vector<vector<double>>>();
                if(dots.size()!=difforder.size()||dots.size()!=value.size()){
                    cout<<"The number of left endpoint sequences, derivative orders, right derivations should be the same."<<endl;
                    throw "invalid shapes";
                }
            }
            if(!value.empty()){
                for (size_t j=0;j<value.size();j++){
                    if((int)value[j].size()<dimension) {cout<<"Lack of complete dimensional information!";throw"Lack of complete dimensional information!";}//检查n阶导数维度
                    values_dimension.push_back(value[j][cur_dimen]);//n阶导值
                }
            }

        }
        if(!j["boundary condition"]["exists"].is_null()){
            json ex =j["boundary condition"]["exists"];
            if(!ex.empty()) exist_dot=ex.get<vector<int>>();
        }
        
        coef.push_back(pp_solve(order,knots,f_value_dimension,dots1,difforder1,dots2,difforder2,dots,difforder,values_dimension,exist_dot));
    }
    

}


vector<double> ppForm::get_value(double x){
    if(x>knots[knots.size()-1]||x<knots[0]||x<begin||x>end){
        cout<<"out of range"<<endl;
        throw "out of range";
    }
    
    int i;
    for (i = (int)knots.size()-2; i >=0; i--)
    {
        if(x>=knots[i]) break;
    }
    vector<double> result;
    for (int k = 0; k < dimension; k++)
    {
        double power=1,sum=0;
        for (int j = 0; j <=order; j++)
        {
            sum+=power*coef[k][i][j];
            power*=x-knots[i];
        }
        result.push_back(sum);
    }
    return result;
}

ppForm::~ppForm()
{
    
}

#endif