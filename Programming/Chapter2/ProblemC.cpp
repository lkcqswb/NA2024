#include<iostream>
#include"PA_Newton_Formula.hpp"
#include<cmath>
#include<vector>
#include <fstream> 
using namespace std;

double f(double x){
    return 1/(1+25*pow(x,2));
}


int main(){
    int n,i,j;
    vector<double> dots;
    ofstream outfile("ProblemC.txt", ios::trunc);
    if (!outfile.is_open()) {
        cerr << "can't open ProblemC.txt" << endl;
        return 1;
    }

    for(j=1;j<5;j++){
        n=5*j;
        for(i=1;i<n+1;i++){
            dots.push_back(cos((2*i-1)*pi/(2*n)));
        }
        Newton_formula New1(f,dots);
        for (i = 0; i < 101; i++) {
            outfile << 0.02*i-1 << " " << New1.get_interpolation_value(0.02*i-1) << endl;
        }
        outfile << "#END# n=" <<n<< endl;
        dots.clear();
    }
    for (i = 0; i < 101; i++) {
            outfile << 0.02*i-1 << " " << f(0.02*i-1) << endl;
    }
    outfile << "#END#f itself"<< endl;


    string command = "python3 plot.py ProblemC.txt";
    system(command.c_str());
    return 0;
}