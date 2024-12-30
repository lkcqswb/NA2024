#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<fstream>
using namespace std;
double truncate(double x,int order){
    if(x>0) return pow(x,order);
    return 0;
}

void test_truncate(double t, string file_name,double begin,double end){
    ofstream outfile(file_name, ios::trunc);
    for (double x = begin; x <= end; x+=(end-begin)/100) outfile << x << " " << truncate(t-x,1) << endl;
    outfile << "#END# " <<"n=1"<< endl;
    for (double x = begin; x <= end; x+=(end-begin)/100) outfile << x << " " << truncate(t-x,2) << endl;
    outfile << "#END# " <<"n=2"<< endl;
}

void divide_truncate(vector<double> t,string file_name,int order){
    if(t.size()<order+2){
        cerr<<"invalid knot"<<endl;
        throw "invalid knot";
    }
    ofstream outfile(file_name, ios::trunc);
    int i=0;
    for (double x = t[0]-1; x <= t[t.size()-1]; x+=0.1)
    {
        vector<double> value(order+2);
        for (i = 0; i < order+2; i++)
            {
                value[i]=truncate(t[i]-x,order);
                outfile<<"t_{"<<i-1<<","<<0<<"}"<< x << " " << value[i] << endl;
            }
        for (size_t layer = 1; layer < order+2; layer++)
        {
            for (i = order+1; i >= layer; i--)
            {
                value[i]=(value[i]-value[i-1])/(t[i]-t[i-layer]);
                outfile<<"t_{"<<i-1<<","<<layer<<"}"<< x << " " << value[i] << endl;
            }
        }
        
    }
}

void plotB_n(vector<double> t,string file_name,int order){
    if(t.size()<order+2){
        cerr<<"invalid knot"<<endl;
        throw "invalid knot";
    }
    ofstream outfile(file_name, ios::trunc);
    int i=0;
    for (double x = t[0]-1; x <= t[t.size()-1]; x+=0.1)
    {
        vector<double> value(order+2);
        for (i = 0; i < order+2; i++)
            {
                value[i]=truncate(t[i]-x,order);
            }
        for (size_t layer = 1; layer < order+2; layer++)
        {
            for (i = order+1; i >= layer; i--)
            {
                value[i]=(value[i]-value[i-1])/(t[i]-t[i-layer]);
                if(layer == order+1) outfile<<x << " " << (t[t.size()-1]-t[0])*value[i] << endl;
            }
        }
    }
    outfile<<"#END#B_0^3"<<endl;
}

int main(){
    test_truncate(0,"PF_fig1.txt",-1,2);
    string command = "python plot.py PF_fig1.txt";
    system(command.c_str());
    divide_truncate({-1,0,1,2},"PF_fig2.txt",1);
    command="python plot_staircase.py PF_fig2.txt";
    system(command.c_str());

    divide_truncate({-1,0,1,2},"PF_fig3.txt",2);
    command="python plot_staircase.py PF_fig3.txt";
    system(command.c_str());

    plotB_n({-1,0,1,2},"PF_fig4.txt",2);
    command="python plot.py PF_fig4.txt";
    system(command.c_str());
}   
