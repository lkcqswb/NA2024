#pragma once
#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/BSpline.hpp"
#include"../../include/Eigen/Dense"
#include"functions.hpp"
using json = nlohmann::json;
using namespace std;



class sphere_fitting_B{
private:
    double start,end,radius;
    int N;
    vector<vector<double>> points;
    BSpline* B_lr;
    BSpline* B_ud;
    vector<double> O;
    
        
public:
    sphere_fitting_B(json j){
        if(!j["range"]["begin"].is_null()) start=j["range"]["begin"];
        else start=0;
        if(!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=1;
        if(start==end){
            cout<<"invalid assignments"<<endl;
            start=0;
            end=1;
        }
        O=j["centre"].get<vector<double>>();
        radius=j["radius"];
        N=j["points"].get<vector<vector<double>>>().size();
       
        vector<double> knots={};
        double r=j["radius"];
        vector<double> centre=j["centre"].get<vector<double>>();
        if(centre.size()<3){
            cout<<"invalid input"<<endl;
            throw "invalid input";
        }
        vector<vector<double>> points_ud,points_lr;

        if(!j["points"].is_null()){
            points=j["points"].get<vector<vector<double>>>();
            points_ud=j["points"].get<vector<vector<double>>>();
            points_lr=j["points"].get<vector<vector<double>>>();
        
            process_value_ud(points_ud,centre,r);
            process_value_lr(points_lr,centre,r);
         
        }else{
            throw "enadequate information";
        }
        for (int i = 0; i < N; ++i) {
            knots.push_back(i);
        }
        vector<int> dots1={},dots2={},difforder1={},difforder2={};
        vector<vector<double>> value1={},value2={};
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
                value1=va[2].get<vector<vector<double>>>();
                value2=va[2].get<vector<vector<double>>>();
                if(dots.size()!=orders.size()||dots.size()!=value1.size()){
                    cout<<"The number of left endpoint sequences, derivative orders, right derivations should be the same."<<endl;
                    throw "invalid shapes";
                }
                for (size_t i = 0; i < value1.size(); i++)
                {
                    value1[i]=process_derivation_ud(points_ud[dots[i]],value1[i],orders[i],r,(double)(N-1)/(end-start));
                    value2[i]=process_derivation_lr(points_lr[dots[i]],value2[i],orders[i],r,(double)(N-1)/(end-start));
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
        value1.reserve(value1.size() + points_ud.size());
        value2.reserve(value2.size() + points_lr.size());
        for (size_t i = 0; i < points.size(); i++) value1.insert(value1.begin(), points_ud[i]); 
        for (size_t i = 0; i < points.size(); i++) value2.insert(value2.begin(), points_lr[i]); 
        json jud={
            {"dimension",2},
            {"order", j["order"]},
            {"boundary condition", {
                {"equals",{dots1,difforder1,dots2,difforder2}},
                {"values", {
                    dots,
                    orders,
                    value1
                }},
                {"exists",exist_dot}
            }},
            {"data points", knots},
            {"range", {
                {"end", N-1},
                {"begin", 0}
            }}
        };
        json jlr={
            {"dimension",2},
            {"order", j["order"]},
            {"boundary condition", {
                {"equals",{dots1,difforder1,dots2,difforder2}},
                {"values", {
                    dots,
                    orders,
                    value2
                }},
                {"exists",exist_dot}
            }},
            {"data points", knots},
            {"range", {
                {"end", N-1},
                {"begin", 0}
            }}
        };
 
        B_lr=new BSpline(jlr);
        B_ud=new BSpline(jud);
    
        
    };

    ~sphere_fitting_B(){
        delete B_lr;
        delete B_ud;
    }

    vector<double> get_value(double t){
        if(t>=start&& t<=end){
            if(abs(points[(int)(t*(N-1))][2])==radius||((int)(t*(N-1))>0&&abs(points[(int)(t*(N-1))-1][2])==radius)||((int)(t*(N-1))<N-1&&abs(points[(int)(t*(N-1))+1][2])==radius)){
                vector<double> angle=B_lr->get_value(t*(N-1));
                double y=radius*sin(angle[0]);
                double refle=radius*cos(angle[0]);
                return {refle*cos(angle[1])+O[0],y+O[1],refle*sin(angle[1])+O[2]};
            }
            vector<double> angle=B_ud->get_value(t*(N-1));
            double z=radius*sin(angle[0]);
            double refle=radius*cos(angle[0]);
            return {refle*cos(angle[1])+O[0],refle*sin(angle[1])+O[1],z+O[2]};
        }
        else return O;
    }
};





