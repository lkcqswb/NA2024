#include<iostream>
#include<cmath>
struct point{
    double x,y;
    point(double x,double y):x(x),y(y){}
};
using namespace std;
int main(){
    const int num=6;
    struct point d[num]={point(0,0),point(1,1),point(1,1),point(1,1),point(2,pow(2,7)),point(2,pow(2,7))};
    double f[6][6];
    for (int i = 0; i < num; i++)
    {
        for (int j = i; j < num; j++)
        {
            if(i==0) f[j][i]=d[j].y;
            else if((d[j].x==d[j-i].x)){
                f[j][i]=(i==1) ? 7*pow(d[j].x,6):21*pow(d[j].x,5);
            }
            else{
                f[j][i]=(f[j][i-1]-f[j-1][i-1])/(d[j].x-d[j-i].x);
            }
        }
        
    }
    for (int j = 0; j < num; j++)
    {
        for (int i = 0; i <= j; i++)
        {
            cout<<f[j][i]<<" ";
        }
        cout<<endl;
    }

    
}