#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/BSpline.hpp"
using json = nlohmann::json;
using namespace std;

class BSpline_1_0: public BSpline{
private:
    vector<double> knots; 
    double start,end;
public:
    BSpline_1_0(json j):BSpline(j){
        if (!j["data points"].is_null()){
                knots=j["data points"].get<vector<double>>();
            if(knots.empty()||knots.size()<2) {cout<<"enadequate knots"<<endl; throw "enadequate datapoints";}
        }else {cout<<"no datapoints"<<endl; throw "no datapoints";}
        
        if (!j["range"]["begin"].is_null()) start=j["range"]["begin"];
        else start=knots[0];
        if (!j["range"]["end"].is_null()) end=j["range"]["end"];
        else end=knots[knots.size()-1];
    };
    double get_value(double t){
        return BSpline::get_value(t)[0];
    }
};





