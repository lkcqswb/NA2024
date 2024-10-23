#include<iostream>
#include"PA_Newton_Formula.hpp"
#include<cmath>
#include <fstream> 
using namespace std;


double Runge(double x){
    return 1/(1+pow(x,2));
}



int main(){
    int i,n,j;
    vector<double> dots;
    Function f(Runge);
    ofstream outfile("ProblemB.txt", ios::trunc);
    if (!outfile.is_open()) {
        cerr << "can't open ProblemB.txt" << endl;
        return 1;
    }

    for(j=1;j<5;j++){
        n=2*j;
        for (i = 0; i < n+1; i++)
        {
            dots.push_back((double)(10*i)/n-5);
        }
        Newton_formula New1(f,dots);
        for (i = 0; i < 101; i++) {
            outfile << 0.1*i-5 << " " << New1.get_interpolation_value(0.1*i-5) << endl;
        }
        outfile << "#END# n=" <<n<< endl;
        dots.clear();
    }
    for (i = 0; i < 101; i++) {
            outfile << 0.1*i-5 << " " << Runge(0.1*i-5) << endl;
        }
    outfile << "#END#The Runge function itself"<< endl;


    string command = "python plot.py ProblemB.txt";
    system(command.c_str());
    return 0;
}