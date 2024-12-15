#include "../../include/json.hpp"
#include "../../ppForm_a_BSpline/BSpline.hpp"
#include<vector>
#include<string>
using json = nlohmann::json;
using namespace std;

vector<double> get_value(vector<double> knots,vector<vector<double>> coef,int order,int dimension,double x){
    int i;
    
    if(x>knots[knots.size()-1]||x<knots[order-1]){
        cout<<"out of range"<<endl;
        throw "out of range";
    }
    if(x==knots[order-1]){//和ppform统一定义域
        i=order-1;
    }
    else{
        for (i = knots.size()-2; i >=order-1; i--)
        {
            if(x>knots[i]) break;
        }
    }
    vector<vector<double>> B= construct_value_table(knots,order,i+1,x);
    vector<double> result;
    for (int k = 0; k < dimension; k++)
    {
        double sum=0;
        for (int j = 0; j <=i+1; j++)
        {
            sum+=B[order][j]*coef[k][j];
        }
        result.push_back(sum);
    }
    return result;
}


int main(){
    json j;
    std::ifstream file("config.json"); 
    if (!file.is_open()) {
            std::cerr << "无法打开文件" << std::endl;
            return 1;
    }
    j = json::parse(file);
    file.close();
    vector<vector<double>> coef=j["coefficients"].get<vector<vector<double>>>();
    vector<double> knots=j["data points"].get<vector<double>>();
    ofstream outfile1("BSpline.txt", ios::trunc);
    int order=j["order"],dimension=j["dimension"];

    vector<double> add_knots={};
    for (size_t i = order-1; i >0; i--) add_knots.push_back(knots[0]-i);
    knots.insert(knots.begin(), add_knots.begin(), add_knots.end());

    for (double t = j["range"]["begin"]; t <= j["range"]["end"]; t+=0.01)
    {
        vector<double> result=get_value(knots,coef,order,dimension,t);
        outfile1 << t << " " ;
        for(size_t i=0;i<result.size();i++) outfile1 << result[i]<< " " ;
        outfile1<<endl;
    }
    string command = "python plot.py BSpline.txt "+to_string(dimension);
    system(command.c_str());





}