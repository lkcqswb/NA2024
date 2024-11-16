#include "ppform_s_1_0.hpp"
#include "BSpline_s_1_0.hpp"
#include <iostream>
#include <fstream>



int main() {
    ppform_1_0 pp({0,1,2,3},{{10},{20},{32},{0}},0,3);
    BSpline_1_0 bs({0,1,2,3},{{10},{20},{32},{0}},0,3);
    size_t i;
    ofstream outfile1("ppform_s_1_0.txt", ios::trunc);
    for (i = 0; i < 101; i++) {
            outfile1 << 0.03*i << " " << pp.get_value(0.03*i) << endl;
        }
    outfile1 << "#END#"<< endl;
    string command = "python3 plot.py ppform_s_1_0.txt";
    system(command.c_str());


    ofstream outfile2("BSpline_s_1_0.txt", ios::trunc);
    for (i = 1; i < 101; i++) {
            outfile2 << 0.03*i << " " << bs.get_value(0.03*i) << endl;
        }
    outfile2 << "#END#"<< endl;
    command = "python3 plot.py BSpline_s_1_0.txt";
    system(command.c_str());


    return 0;

}