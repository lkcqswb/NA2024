#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/ppForm.hpp"
#include<vector>
using json = nlohmann::json;
using namespace std;

vector<int> generate_seq(int n){
    vector<int> result;
    for (size_t i = 0; i < n; i++)
    {
        result.push_back(i);
    }
    return result;
}
vector<vector<double>> generate_func_value(vector<double> f_values){
    vector<vector<double>> result;
    for (size_t i = 0; i < f_values.size(); i++)
    {
        result.push_back({f_values[i]});
    }
    return result;
}



class ppform_2_3_periodic: public ppForm{
private:
    vector<double> knots; 
    double start,end;
    vector<int> dir,dot;
public:
    ppform_2_3_periodic(vector<double>input_knots,vector<double> f_values,double istart,double iend):
    ppForm({
            {"dimension", 1},
            {"order", 3},
            {"boundary condition", {
                {"equals", {
                    {0, 0, 0},
                    {0, 1, 2},
                    {static_cast<int>(input_knots.size()) - 1, static_cast<int>(input_knots.size()) - 1, static_cast<int>(input_knots.size()) - 1},
                    {0, 1, 2}
                }},
                {"values", {
                    generate_seq(input_knots.size() - 1),
                    vector<int>(input_knots.size() - 1, 0), 
                    generate_func_value(f_values)
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
        return ppForm::get_value(t)[0];
    }
};




