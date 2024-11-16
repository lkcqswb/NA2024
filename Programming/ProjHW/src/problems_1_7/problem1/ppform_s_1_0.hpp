#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/ppForm.hpp"
using json = nlohmann::json;
using namespace std;

class ppform_1_0: public ppForm{
private:
    vector<double> knots;
    vector<vector<double>> f_values;  
    double start,end;
public:
    ppform_1_0(vector<double> i_knots,vector<vector<double>> i_f_values,double astart,double aend):ppForm({
        {"dimension", 1},
        {"order", 1},
        {"boundary condition",{
            {"values",{}},
            {"equals",{}},
            {"exists",{}}
        }},
        {"data points", {i_knots, i_f_values}},
        {"range", {
            {"end", aend},
            {"start", astart}
        }}
    }),knots(i_knots),f_values(i_f_values),start(astart),end(aend){};
    double get_value(double t){
        return ppForm::get_value(t)[0];
    }
};




