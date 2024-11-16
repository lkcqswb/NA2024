#include "ppForm_2_3_periodic.hpp"
#include "ppForm_2_3_complete.hpp"
#include "ppForm_2_3_natural.hpp"
#include<iostream>
#include <vector>
using namespace std;

int main() {
    double start=0,end=10;
    vector<double> knots={0,2,6,7,10};
    vector<double> f_values={4,-7,1,4};
    ppform_2_3_periodic pp(knots,f_values,start,end);
    f_values={4,-7,1,4,7};
    ppform_2_3_natural pn(knots,f_values,start,end);
    vector<double> deri={1,2};
    ppform_2_3_complete pc(knots,f_values,deri,start,end);

    double t;
    ofstream outfile1("ppform_2_3_periodic.txt", ios::trunc);
    ofstream outfile2("ppform_2_3_natural.txt", ios::trunc);
    ofstream outfile3("ppform_2_3_complete.txt", ios::trunc);

    for (t = start; t <= end; t+=0.01) {
            outfile1 << t << " " << pp.get_value(t) << endl;
            outfile2 << t << " " << pn.get_value(t) << endl;
            outfile3 << t << " " << pc.get_value(t) << endl;
    }



    string command = "python plot.py ppform_2_3_periodic.txt";
    system(command.c_str());
    command = "python plot.py ppform_2_3_natural.txt";
    system(command.c_str());
    command = "python plot.py ppform_2_3_complete.txt";
    system(command.c_str());


    return 0;

}