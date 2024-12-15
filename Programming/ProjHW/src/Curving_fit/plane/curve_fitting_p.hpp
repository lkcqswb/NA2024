#pragma once

#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/ppForm.hpp"
using json = nlohmann::json;
using namespace std;


class curve_fitting_p: public ppForm{
private:
    double start,end;
    json construct(json j,string sampling_mode){
        vector<double> knots={};
        int dimension;
        if (!j["dimension"].is_null()) dimension=j["dimension"];
        else cerr<<"no dimension input"<<endl;

        vector<vector<double>> points;
        if(!j["points"].is_null()){
            points=j["points"].get<vector<vector<double>>>();
        }else{
            throw "enadequate information";
        }
        int N=points.size();
        if (!j["range"]["begin"].is_null()) start=j["range"]["begin"];
        else start=0;
        if (!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=1;
        
        if(sampling_mode=="equal"){
            for (int i = 0; i < N; ++i) {
                double t = static_cast<double>(i) / (N - 1); 
                knots.push_back((end-start)*t+start);
            }
        }else if(sampling_mode=="Chord"){
            vector<double> chord={0};
            for (size_t i = 1; i <points.size(); i++)//计算累计弦长
            {
                double len=0;
                for (size_t j = 0; j < points[i].size(); j++) len+=(points[i][j]-points[i-1][j])*(points[i][j]-points[i-1][j]);
                len=sqrt(len)+chord[chord.size()-1];
                chord.push_back(len);
            }
            for (size_t i = 0; i < points.size(); i++)
            {
                knots.push_back((end-start)*chord[i]/chord[chord.size()-1]+start);
            }
        }else{
            cerr<<"no such sampling_mode"<<endl;
            throw "no such sampling_mode";
        }
        
        
        vector<int> dots1={},dots2={},difforder1={},difforder2={};
        vector<vector<double>> value={};
        vector<int> dots={},orders={};
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
            if (!va.empty()){
                if(va.size()<3){
                    cout<<"illegal boundary condition"<<endl;
                    throw "illegal boundary condition";
                }
                dots=va[0].get<vector<int>>();
                orders=va[1].get<vector<int>>();
                value=va[2].get<vector<vector<double>>>();
                if(dots.size()!=orders.size()||dots.size()!=value.size()){
                    cout<<"The number of left endpoint sequences, derivative orders, right derivations should be the same."<<endl;
                    throw "invalid shapes";
                }
            }
        }
        if(!j["boundary condition"]["exists"].is_null()){
            json ex =j["boundary condition"]["exists"];
            if(!ex.empty()) exist_dot=ex.get<vector<int>>();
        }
    
        dots.reserve(dots.size() + knots.size());
        for (size_t i = 0; i < knots.size(); i++) dots.insert(dots.begin(), i); 
        

        // 预留空间为 orders 添加零以避免重新分配
        orders.reserve(orders.size() + knots.size());
        orders.insert(orders.begin(), knots.size(), 0); // 在 orders 的末尾插入 knots.size() 个零

        // 预留空间以避免重新分配
        value.reserve(value.size() + points.size());
        value.insert(value.begin(), points.begin(), points.end());


        return {
            {"dimension",dimension},
            {"order", j["order"]},
            {"boundary condition", {
                {"equals",{dots1,difforder1,dots2,difforder2}},
                {"values", {
                    dots,
                    orders,
                    value
                }},
                {"exists",exist_dot}
            }},
            {"data points", knots},
            {"range", {
                {"end", end},
                {"begin", start}
            }}
        };
        
    }
public:
    curve_fitting_p(json j,string sampling_mode="equal"):ppForm(construct(j,sampling_mode)){
        if(!j["range"]["begin"].is_null()) start=j["range"]["begin"];
        else start=0;
        if(!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=1;

    };
    vector<double> get_value(double t){
        if(t>=start&& t<=end) return ppForm::get_value(t);
        else return {0,0};
    }
};





