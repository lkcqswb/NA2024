#include"Newton_Formula_without_F.hpp"
#include<iostream>
#include<vector>
#include <fstream> 
using namespace std;


int main(){
    vector<vector<double>> dots1={{0,6.67},{6,17.3},{10,42.7},{13,37.3},{17,30.1},{20,29.3},{28,28.7}};
    vector<vector<double>> dots2={{0,6.67},{6,16.1},{10,18.9},{13,15.0},{17,10.6},{20,9.44},{28,8.89}};
    Newton_formula_WF  New1(dots1);
    Newton_formula_WF  New2(dots2);
    ofstream outfile("ProblemE.txt", ios::trunc);
    int day;
    if (!outfile.is_open()) {
        cerr << "can't open ProblemE.txt" << endl;
        return 1;
    }
    for (day = 0; day < 30; day++) {
            outfile << day << " " << New1.get_interpolation_value(day) << endl;
    }
    outfile << "#END# Sp1" << endl;
     for (day = 0; day < 30; day++) {
            outfile << day << " " << New2.get_interpolation_value(day) << endl;
    }
    outfile << "#END# Sp2" << endl;   
    cout<<"x represents time, and y represents weight."<<endl;
    string command = "python plot.py ProblemE.txt";
    system(command.c_str());

    cout << "Sp1 after another 15 days: average weight: " << New1.get_interpolation_value(43)  << endl;
    cout<<"Abnormal weight, predicted to die."<<endl;
    cout << "Sp2 after another 15 days: average weight: " << New2.get_interpolation_value(43)  << endl;
    cout<<"Abnormal weight, predicted to die."<<endl;
}