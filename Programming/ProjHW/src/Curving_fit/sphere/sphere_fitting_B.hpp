#pragma once
#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/BSpline.hpp"
#include"../../include/Eigen/Dense"
#include"functions.hpp"
using json = nlohmann::json;
using namespace std;



class sphere_fitting_B: public BSpline{
private:
    double start,end,radius;
    int N;
    vector<double> O;
    json construct(json j){
        vector<double> knots={};
        double r=j["radius"];
        vector<double> centre=j["centre"].get<vector<double>>();
        if(centre.size()<3){
            cout<<"invalid input"<<endl;
            throw "invalid input";
        }
        vector<vector<double>> points;
        if(!j["points"].is_null()){
            points=j["points"].get<vector<vector<double>>>();
            for (size_t i = 0; i < points.size(); i++) points[i]=process_value(points[i],centre,r);
        }else{
            throw "enadequate information";
        }
        int N=points.size();

        for (int i = 0; i < N; ++i) {
            knots.push_back(i);
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
                for (size_t i = 0; i < value.size(); i++)
                {
                    value[i]=process_derivation(points[dots[i]],value[i],orders[i],r);
                }
                
            }
        }
        if(!j["boundary condition"]["exists"].is_null()){
            json ex =j["boundary condition"]["exists"];
            if(!ex.empty()) exist_dot=ex.get<vector<int>>();
        }
    
        dots.reserve(dots.size() + knots.size());
        for (size_t i = 0; i < knots.size(); i++) dots.insert(dots.begin(), i); 
        orders.reserve(orders.size() + knots.size());
        orders.insert(orders.begin(), knots.size(), 0); 
        value.reserve(value.size() + points.size());
        for (size_t i = 0; i < points.size(); i++) value.insert(value.begin(), points[i]); 

        return {
            {"dimension",2},
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
                {"end", N-1},
                {"begin", 0}
            }}
        };
        
    }
public:
    sphere_fitting_B(json j):BSpline(construct(j)){
        if(!j["range"]["begin"].is_null()) start=j["range"]["begin"];
        else start=0;
        if(!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=1;
        O=j["centre"].get<vector<double>>();
        radius=j["radius"];
        N=j["points"].get<vector<vector<double>>>().size();
    };
    vector<double> get_value(double t){
        if(t>=start&& t<=end){
            vector<double> angle=BSpline::get_value(t*(N-1));
            double z=radius*sin(angle[0]);
            double refle=radius*cos(angle[0]);
            return {refle*cos(angle[1])+O[0],refle*sin(angle[1])+O[1],z+O[2]};
        }
        else return O;
    }
};





