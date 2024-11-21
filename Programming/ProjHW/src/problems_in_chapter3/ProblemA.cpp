#include "..\problems_1_7\problem2\ppform_2_3_natural.hpp"
#include "..\problems_1_7\problem3\BSpline_2_3_natural.hpp"
#include<iostream>
#include<vector>
using namespace std;
vector<double> construct(int N,double start,double end){
    vector<double> knots;
    for (int i = 0; i < N; i++)
    {
        knots.push_back((double)i/(N-1)*(end-start)+start);
    }
    return knots;
}

double f(double x){
    return 1/(1+25*x*x);
}

vector<double> construct_f(vector<double> knots){
    vector<double> values;
    for (size_t i = 0; i < knots.size(); i++)
    {
        values.push_back(f(knots[i]));
    }
    return values;
}


double test_max_norm_p(int N,double start,double end){
    vector<double> knots=construct(N,start,end);
    vector<double> f_values=construct_f(knots);
    ppform_2_3_natural pp(knots,f_values,start,end);
    double max_error=0;
    for (size_t i = 0; i < knots.size()-1; i++)
    {
        double result=pp.get_value((knots[i]+knots[i+1])/2);
        double result2=f((knots[i]+knots[i+1])/2);
        if(abs(result2-result)>max_error) max_error=abs(result2-result);
    }
    return max_error;
}

double test_max_norm_B(int N,double start,double end){
    vector<double> knots=construct(N,start,end);
    vector<double> f_values=construct_f(knots);
    BSpline_2_3_natural bs(knots,f_values,start,end);
    double max_error=0;
    for (size_t i = 0; i < knots.size(); i++)
    {
        double result=bs.get_value((knots[i]+knots[i+1])/2);
        double result2=f((knots[i]+knots[i+1])/2);
        if(abs(result2-result)>max_error) max_error=abs(result2-result);
    }
    return max_error;
}

int main(){
    int n[5]={6,11,21,41,81};
    for (size_t i = 0; i < 5; i++)
    {
        cout<<n[i]<<": "<<endl<<test_max_norm_p(n[i],-1,1)<<endl;
        cout<<n[i]<<": "<<endl<<test_max_norm_B(n[i],-1,1)<<endl;
    }
    

}

