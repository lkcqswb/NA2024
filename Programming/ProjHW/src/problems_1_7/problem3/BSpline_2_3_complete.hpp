#pragma once
#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/BSpline.hpp"
#include<vector>
using json = nlohmann::json;
using namespace std;


class BSpline_2_3_complete: public BSpline{
private:

    double start,end;

    json process_json(json j){
        if (!j.contains("data points") || !j["data points"].is_array()) {
            cerr << "Error: 'data points' is missing" << endl;
            throw "Invalid";
        }

        vector<double> knots;
        for (const auto& item : j["data points"]) {
            if (item.is_number()) {
                knots.push_back(item.get<double>());
            } else {
                cerr << "Error: Non-numeric value in 'data points'" << endl;
                throw "Non-numeric value in 'data points'";
            }
        }

        if (knots.size() < 2) {
            cerr << "Error: 'data points' must contain at least two values" << endl;
            throw "Insufficient data points";
        }

        vector<double> values;
        if (j.contains("function values") && j["function values"].is_array()) {
            for (const auto& item : j["function values"]) {
                if (item.is_number()) {
                    values.push_back(item.get<double>());
                } else {
                    cerr << "Error: Non-numeric value in 'function values'" << endl;
                    throw "Non-numeric value in 'function values'";
                }
            }
        } else {
            cerr << "Error: 'function values' is missing or not an array" << endl;
            throw "'function values' missing or invalid";
        }
        if(values.size()!=knots.size()){
            throw "invalid";
        }

        if (!j["range"]["begin"].is_null()) start=j["range"]["begin"];
        else start=knots[0];
        if (!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=knots[knots.size()-1];

        vector<double> value;
        if (!j["derivation"].is_null()) value=j["derivation"].get<vector<double>>();
        if(j["derivation"].is_null()||value.size()!=2){
            throw "invalid";
        }


        vector<double> points,order;
        points.push_back(0);
        order.push_back(1);
        for (size_t i = 0; i < knots.size(); i++)
        {
            points.push_back(i);
            order.push_back(0);
        }

        points.push_back(knots.size()-1);
        order.push_back(1);
        values.insert(values.begin(), value[0]);
        values.push_back(value[1]);

        vector<vector<double>> input_value;
        for (size_t i = 0; i < values.size(); i++)
        {
            input_value.push_back({values[i]});
        }
        return json {
            {"dimension", 1},
            {"order", 3},
            {"boundary condition", {
                {"values", {
                    points,
                    order,
                    input_value
                }},
            }},
            {"data points", knots},
            {"range", {
                {"end", end},
                {"begin", start}
            }}
        };
        
    }
        
public:
    BSpline_2_3_complete(json j):BSpline(process_json(j)){};
    double get_value(double t){
        return BSpline::get_value(t)[0];
    }
};







//因为题目C中要求使用Theorem 3.57/3.58.而上面这个思路不同，于是决定重新实现一遍
class BSpline_uniform_23{
private:
    int max(int x,int y){
        if(x>y) return x;
        return y;
    }

    vector<double> update_B(vector<double> B,vector<double> knots,int current_order,double x){
        for (int k = max(0,B.size()-2-current_order); k <= (int)B.size()-2; k++)
        {
            B[k]=B[k]*(x-get_ti(k-1,knots))/(get_ti(current_order+k,knots)-get_ti(k-1,knots))+B[k+1]*(get_ti(k+current_order+1,knots)-x)/(get_ti(current_order+k+1,knots)-get_ti(k,knots));
        }
        B[B.size()-1]*=(x-get_ti(B.size()-2,knots))/(get_ti(current_order+B.size()-1,knots)-get_ti(B.size()-2,knots));
        return B;
    }

    vector<vector<double>> construct_value_table(vector<double> knots,int order,int index,double x){
        vector<vector<double>> table={};
        vector<double> B(index+1,0);
        B[index]=1;
        table.push_back(B);
        for (int i = 0; i < order; i++)
        {
            table.push_back(update_B(table[table.size()-1],knots,i,x));
        }
        return table;
    }
    double get_ti(int index,vector<double> knots){
        double intervals=(knots[knots.size()-1]-knots[0])/knots.size();
        if(index<0) return knots[0]+index*intervals;
        if(index<=(int)knots.size()-1) return knots[index];
        else return knots[knots.size()-1]+(-knots.size()+1+index)*intervals;
    }

    vector<double> knots,new_knots;
    vector<double> coeff;
    double begin,end,interval;
public:
    BSpline_uniform_23(json j){
        
        knots=j["data points"].get<vector<double>>();
        if (!j["range"]["begin"].is_null()) begin=j["range"]["begin"];
        else begin=knots[0];
        if (!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=knots[knots.size()-1];


        vector<double> value=j["values"].get<vector<double>>();
         vector<double> diff_value=j["diff_values"].get<vector<double>>();
        Eigen::MatrixXd m=Eigen::MatrixXd::Zero(knots.size(),knots.size());
        interval = (knots[knots.size()-1]-knots[0])/(knots.size()-1);
        new_knots={};
        for (size_t i = 0; i < knots.size(); i++)
        {
            new_knots.push_back(knots[i]/interval);
        }
        new_knots.insert(new_knots.begin(),new_knots[0]-1);
        new_knots.insert(new_knots.begin(),new_knots[0]-1);

        for (size_t i = 0; i < knots.size(); i++)
        {
            if(i==0){
                m(0,0)=2;
                m(0,1)=1;
            }else if(i==knots.size()-1){
                m(knots.size()-1,knots.size()-2)=1;
                m(knots.size()-1,knots.size()-1)=2;
            }else{
                m(i,i)=4;
                m(i,i-1)=1;
                m(i,i+1)=1;
            }
        }
        Eigen::VectorXd target(knots.size());
        for (size_t i = 0; i < knots.size(); i++)
        {
            if(i==0){
                target[i]=3*value[i]+diff_value[0]*interval;
            }else if(i==knots.size()-1){
                target[i]=3*value[knots.size()-1]-diff_value[1]*interval;
            }else{
                target[i]=6*value[i];
            }
        }
        Eigen::VectorXd solution = m.colPivHouseholderQr().solve(target);

        coeff={solution[1]-2*diff_value[0]*interval};

        for (size_t i = 0; i < solution.size(); i++)
        {
            coeff.push_back(solution[i]);
        }

        coeff.push_back(solution[solution.size()-2]+2*diff_value[1]*interval);
        knots.insert(knots.begin(),knots[0]-interval);
        knots.insert(knots.begin(),knots[0]-interval);


    }
    double get_value(double x){
        int i;
        if(x>knots[knots.size()-1]||x<knots[2]||x<begin||x>end){
            cout<<x<<" is out of range"<<endl;
            throw "out of range";
        }
        if(x==knots[2]){//和ppform统一定义域
            i=2;
        }
        else{
            for (i = knots.size()-2; i >=0; i--)
            {
                if(x>knots[i]) break;
            }
        }
        x/=interval;
        vector<vector<double>> B= construct_value_table(new_knots,3,i+1,x);
        double sum=0;
            for (int j = 0; j <=i+1; j++)
            {
                sum+=B[3][j]*coeff[j];
            }
        return sum;
    }
};




class BSpline_uniform_12{
private:
    int max(int x,int y){
        if(x>y) return x;
        return y;
    }

    vector<double> update_B(vector<double> B,vector<double> knots,int current_order,double x){
        for (int k = max(0,B.size()-2-current_order); k <= (int)B.size()-2; k++)
        {
            B[k]=B[k]*(x-get_ti(k-1,knots))/(get_ti(current_order+k,knots)-get_ti(k-1,knots))+B[k+1]*(get_ti(k+current_order+1,knots)-x)/(get_ti(current_order+k+1,knots)-get_ti(k,knots));
        }
        B[B.size()-1]*=(x-get_ti(B.size()-2,knots))/(get_ti(current_order+B.size()-1,knots)-get_ti(B.size()-2,knots));
        return B;
    }

    vector<vector<double>> construct_value_table(vector<double> knots,int order,int index,double x){
        vector<vector<double>> table={};
        vector<double> B(index+1,0);
        B[index]=1;
        table.push_back(B);
        for (int i = 0; i < order; i++)
        {
            table.push_back(update_B(table[table.size()-1],knots,i,x));
        }
        return table;
    }

    double get_ti(int index,vector<double> knots){
        double interval=(knots[knots.size()-1]-knots[0])/knots.size();
        if(index<0) return knots[0]+index*interval;
        if(index<=(int)knots.size()-1) return knots[index];
        else return knots[knots.size()-1]+(-knots.size()+1+index)*interval;
    }
    vector<double> knots,new_knots;
    vector<double> coeff;
    double begin,end,interval;

public:
    BSpline_uniform_12(json j){
        knots=j["data points"].get<vector<double>>();//N-1个元素
        vector<double> value=j["values"].get<vector<double>>();
        vector<double> boundary_values=j["boundary_values"].get<vector<double>>();

        Eigen::MatrixXd m=Eigen::MatrixXd::Zero(knots.size(),knots.size());
        interval = (knots[knots.size()-1]-knots[0])/(knots.size()-1);

        for (size_t i = 0; i < knots.size(); i++)
        {
            if(i==0){
                m(0,0)=5;
                m(0,1)=1;
            }else if(i==knots.size()-1){
                m(knots.size()-1,knots.size()-1)=5;
                m(knots.size()-1,knots.size()-2)=1;
            }else{
                m(i,i)=6;
                m(i,i-1)=1;
                m(i,i+1)=1;
            }
        }
        Eigen::VectorXd target(knots.size());
        for (size_t i = 0; i < knots.size(); i++)
        {
            if(i==0){
                target[i]=8*value[i]-2*boundary_values[0];
            }else if(i==knots.size()-1){
                target[i]=8*value[knots.size()-1]-2*boundary_values[1];
            }else{
                target[i]=8*value[i];
            }
        }
        Eigen::VectorXd solution = m.colPivHouseholderQr().solve(target);
        
        coeff={-solution[0]+2*boundary_values[0]};

        for (size_t i = 0; i < solution.size(); i++)
        {
            coeff.push_back(solution[i]);
        }
        coeff.push_back(-solution[solution.size()-1]+2*boundary_values[1]);

        if (!j["range"]["begin"].is_null()) begin=j["range"]["begin"];
        else begin=knots[0];
        if (!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=knots[knots.size()-1];
         //增加结点

        double interval = (knots[knots.size()-1]-knots[0])/(knots.size()-1);
        knots.insert(knots.begin(), knots[0]-interval);
        knots.push_back(knots[knots.size()-1]+interval);
        for (size_t i = 0; i < knots.size(); i++)
        {
            new_knots.push_back(knots[i]/interval);
        }

    }
    double get_value(double x){
        int i;
        if(x>knots[knots.size()-2]||x<knots[1]||x<begin||x>end){
            cout<<x<<" is out of range"<<endl;
            throw "out of range";
        }
        if(x==knots[1]){
            i=1;
        }
        else{
            for (i = knots.size()-2; i >=1; i--)
            {
                if(x>knots[i]) break;
            }
        }
        x/=interval;
        vector<vector<double>> B= construct_value_table(new_knots,2,i+1,x);
        double sum=0;
            for (int j = 0; j <=i+1; j++)
            {
                sum+=B[2][j]*coeff[j];
            }
        return sum;
    }
};