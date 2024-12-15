#pragma once
#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/ppForm.hpp"
#include<vector>
using json = nlohmann::json;
using namespace std;


class ppform_2_3_complete: public ppForm{
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
    ppform_2_3_complete(json j):ppForm(process_json(j)){};
    double get_value(double t){
        return ppForm::get_value(t)[0];
    }
};




