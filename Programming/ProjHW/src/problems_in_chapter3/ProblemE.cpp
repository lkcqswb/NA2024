#include "../Curving_fit/plane/curve_fitting_B.hpp"
#include "../Curving_fit/plane/curve_fitting_p.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "../include/json.hpp"
#include "../bonus/intersect.hpp"
#define pi 3.1415926535
using json = nlohmann::json;


//t 为0-1
vector<double>r1(double t){
    double x_temp=sqrt(3)*cos(t*2*pi);
    double y_temp=2*(sqrt(sqrt(3)*abs(cos(2*pi*t)))-sqrt(3)*sin(2*pi*t))/3;
    return {x_temp,y_temp};
}
//t 为0-6pi
vector<double>r2(double t){
    double x_temp=sin(t)+t*cos(t);
    double y_temp=cos(t)-t*sin(t);
    return {x_temp,y_temp};
}

//t为0-2pi
vector<double>r3(double t){
    double x_temp=sin(cos(t))*cos(sin(t));
    double y_temp=sin(cos(t))*sin(sin(t));
    double z_temp=cos(cos(t));
    return {x_temp,y_temp,z_temp};
}




void write_errors_to_file(int num_control,int num_plot) {
    ofstream outfile("Pe.txt");
    

    // 控制点的生成,题目中没规定控制点的选取，就随便选了，参数化在后续程序中
    vector<vector<double>> knot1, knot2, knot3;
    for (size_t i = 0; i <= num_control; i++) {
        knot1.push_back(r1((double)(i) / num_control));
        knot2.push_back(r2((double)(i * 6 * pi) / num_control));
        knot3.push_back(r3((double)(i * 2 * pi) / num_control));
    }

    // JSON 格式的参数设置
    json j1 = {
        {"dimension", 2},
        {"order", 3},
        {"boundary condition", {{"exists", {1, knot1.size() - 2}}}},
        {"points", knot1},
        {"range", {{"end", 1}, {"begin", 0}}}
    };
    json j2 = {
        {"dimension", 2},
        {"order", 3},
        {"boundary condition", {{"exists", {1, knot2.size() - 2}}}},
        {"points", knot2},
        {"range", {{"end", 6*pi}, {"begin", 0}}}
    };
    json j3 = {
        {"dimension", 3},
        {"order", 3},
        {"boundary condition", {{"exists", {1, knot3.size() - 2}}}},
        {"points", knot3},
        {"range", {{"end", 2*pi}, {"begin", 0}}}
    };

    cout<<"fig1:intersect?"<<": "<<detect_intersect(j1)<<endl;
    cout<<"fig2:intersect?"<<": "<<detect_intersect(j2)<<endl;
    cout<<"fig3:intersect?"<<": "<<detect_intersect(j3)<<endl;
    // 计算并创建 curve_fitting 对象
    string str1 = "equal", str2 = "Chord";
    curve_fitting_p pe1(j1, str1), pe2(j2, str1), pe3(j3, str1);
    curve_fitting_B be1(j1, str1), be2(j2, str1), be3(j3, str1);
    curve_fitting_p pc1(j1, str2), pc2(j2, str2), pc3(j3, str2);
    curve_fitting_B bc1(j1, str2), bc2(j2, str2), bc3(j3, str2);

    // 计算并写入文件
    cout<<"fig1 "<<endl;
    outfile << "fig1 "<<endl;

    for (size_t i = 0; i <= num_plot; i++) {
        double t1 = min((double)(i) / num_plot,1.0);
        vector<double> pe1_value = pe1.get_value(t1);
        vector<double> be1_value = be1.get_value(t1);
        vector<double> pc1_value = pc1.get_value(t1);
        vector<double> bc1_value = bc1.get_value(t1);

        outfile <<" equal length ppForm,"<<pe1_value[0] << "," << pe1_value[1] <<endl;
        outfile <<" equal length BSpline,"<<be1_value[0] << "," << be1_value[1] << endl;

        outfile << "cumulative chordal length ppform," << pc1_value[0] << "," << pc1_value[1] << endl;
        outfile << "cumulative chordal length BSpline," << bc1_value[0] << "," << bc1_value[1] << endl;
    }
    cout<<"fig2 "<<endl;
    outfile << "fig2 "<<endl;
    for (size_t i = 0; i <= num_plot; i++) {
        double t2 = min((double)(6 * pi * i) / num_plot,6*3.1415);
        vector<double> pe2_value = pe2.get_value(t2);
        vector<double> be2_value = be2.get_value(t2);
        vector<double> pc2_value = pc2.get_value(t2);
        vector<double> bc2_value = bc2.get_value(t2);

        outfile <<" equal length ppForm,"<<pe2_value[0] << "," << pe2_value[1] <<endl;
        outfile <<" equal length BSpline,"<<be2_value[0] << "," << be2_value[1] << endl;

        outfile << "cumulative chordal length ppform," << pc2_value[0] << "," << pc2_value[1] << endl;
        outfile << "cumulative chordal length BSpline," << bc2_value[0] << "," << bc2_value[1] << endl;
    }
    cout<<"fig3 "<<endl;
    outfile << "fig3 "<<endl;
    for (size_t i = 0; i <= num_plot; i++) {
        try{
            double t3 = min((double)(2 * pi * i) / num_plot,6.28315);
            vector<double> pe3_value = pe3.get_value(t3);
            vector<double> be3_value = be3.get_value(t3);
            vector<double> pc3_value = pc3.get_value(t3);
            vector<double> bc3_value = bc3.get_value(t3);
            outfile <<" equal length ppForm,"<<pe3_value[0] << "," << pe3_value[1]<<","<<pe3_value[2] <<endl;
            outfile <<" equal length BSpline,"<<be3_value[0] << "," << be3_value[1]<<","<<be3_value[2]<< endl;

            outfile << "cumulative chordal length ppform," << pc3_value[0] << "," << pc3_value[1]<<","<<pc3_value[2] << endl;
            outfile << "cumulative chordal length BSpline," << bc3_value[0] << "," << bc3_value[1]<<","<<bc3_value[2] << endl;
        }
        catch (const char* e) {
        }
    }
    
        
    outfile.close();

}

void call_python_script(string file_name) {
    const char* pythonScript = "plotE.py";
    system(("python " + string(pythonScript)+" "+file_name).c_str());
}

int main() {
    write_errors_to_file(10,1000);  
    call_python_script("pE_10");   
    write_errors_to_file(40,1000);  
    call_python_script("pE_40");  
    write_errors_to_file(160,160);  
    call_python_script("pE_160");  
    return 0;
}

