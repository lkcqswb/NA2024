#ifndef PPFORM
#define PPFORM
#include<iostream>
#include "../include/json.hpp"
#include "ppForm_options/solve_boundary_condition.hpp"
#include "ppForm_options/solve_boundary_condition_version2.hpp"
#include "ppForm_options/solve_boundary_condition_version3.hpp"
#include<fstream>
#include<string>
#include<vector>
#include <algorithm>
using namespace std;
using json = nlohmann::json;
class ppForm
{
private:
    vector<double> knots,new_knots;
    double interval;
    vector<vector<vector<double>>> coef;
    double begin,end;
    int order;
    int dimension;
public:
    int version;
    ppForm(json j,int version=1);
    ~ppForm();
    vector<double> get_value(double x);
};


ppForm::ppForm(json j,int version)
{
    
    //检查输入的json是否合法
    if (!j["dimension"].is_null()) dimension=j["dimension"];
    if (!j["order"].is_null()) order=j["order"];
    if (!j["data points"].is_null()){
        knots=j["data points"].get<vector<double>>();
        if(knots.empty()||knots.size()<2) {cout<<"enadequate knots"<<endl; throw "enadequate datapoints";}
    }else {cout<<"no datapoints"<<endl; throw "no datapoints";}
    
    //对结点排序
    sort(knots.begin(), knots.end(), [](const double& a, const double& b) {
        return a < b;
        });

    //排序结束



    if (!j["range"]["begin"].is_null()) begin=j["range"]["begin"];
    else begin=knots[0];
    if (!j["range"]["end"].is_null()) end=j["range"]["end"];
    else end=knots[knots.size()-1];
    //结束检查

    //缩放结点
    new_knots={};
    interval = (knots[knots.size()-1]-knots[0])/((knots.size()-1));
    //interval=2;
    for (size_t i = 0; i < knots.size(); i++)
    {
        new_knots.push_back(knots[i]/interval);
    }
    

    
    for(int cur_dimen=0;cur_dimen<dimension;cur_dimen++){


             
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
                    values_dimension.push_back(value[j][cur_dimen]*pow(interval,difforder[j]));//n阶导值
                }
            }

        }
        if(!j["boundary condition"]["exists"].is_null()){
            json ex =j["boundary condition"]["exists"];
            if(!ex.empty()) exist_dot=ex.get<vector<int>>();
        }
        if(version==2){
            coef.push_back(pp_solve_2(order,new_knots,dots1,difforder1,dots2,difforder2,dots,difforder,values_dimension,exist_dot));
        }
        else if(version==3){
            coef.push_back(pp_solve_3(order,new_knots,dots1,difforder1,dots2,difforder2,dots,difforder,values_dimension,exist_dot));
        }else{
            coef.push_back(pp_solve(order,new_knots,dots1,difforder1,dots2,difforder2,dots,difforder,values_dimension,exist_dot));
        }
        
    }


}


vector<double> ppForm::get_value(double x){
    if(x>knots[knots.size()-1]||x<knots[0]||x<begin||x>end){
        if(abs(x-min(end,knots[knots.size()-1]))<1e-4){
            x=min(end,knots[knots.size()-1]);
        }else if(abs(x-max(begin,knots[0]))<1e-4){
            x=max(begin,knots[0]);
        }
        else{
            cout<<x<<" is out of range"<<endl;
            throw std::runtime_error( "out of range");
        }
    }
    
    int i;
    for (i = (int)knots.size()-2; i >=0; i--)
    {
        if(x>=knots[i]) break;
    }
    vector<double> result;
    x/=interval;
    for (int k = 0; k < dimension; k++)
    {
        double power=1,sum=0;
        for (int j = 0; j <=order; j++)
        {
            sum+=power*coef[k][i][j];
            power*=x-new_knots[i];
        }
        result.push_back(sum);
    }
    return result;
}

ppForm::~ppForm()
{
    
}

#endif