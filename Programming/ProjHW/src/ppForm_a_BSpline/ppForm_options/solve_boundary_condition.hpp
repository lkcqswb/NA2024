#include<iostream>
#include<vector>
#include<cmath>
#include"../../include/Eigen/Dense"
#include"ppform_utils.hpp"
using namespace std;



void add_diff(Eigen::MatrixXd& matrix,int line_index,int order,vector<double> input_knots,int dot_index,int difforders,double coeff){
    double delta=input_knots[dot_index]-input_knots[dot_index-1];
    for (int k = (order+1)*(dot_index-1)+difforders; k < (order+1)*(dot_index); k++)
        {
            int offset=k-(order+1)*(dot_index-1);
            matrix(line_index,k)=pow(delta,offset-difforders)*get_fas(offset)/get_fas(offset-difforders)*coeff;
        }
}





vector<vector<double>> pp_solve(int order,vector<double> input_knots,vector<int> dots1,vector<int> difforder1,vector<int> dots2,vector<int> difforder2,vector<int> dots,vector<int> difforder,vector<double> value,vector<int> exist){
    int condition_numbers=dots1.size()+dots.size()+exist.size();//共额外条件数
    int N=input_knots.size();
    int result_size=(N-1)*(order+1);
    if(condition_numbers<N+order-1){
        cout<<"enadequate information"<<endl;
        throw "enadequate information";
    }
    

    //构建方程
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(condition_numbers+order*N-2*order, result_size);
    Eigen::VectorXd target(condition_numbers+order*N-2*order);
    target.setZero();
    size_t line_index=order*N-2*order-1;
    for (int t=1; t < N-1; t++)
    {
        double delta=input_knots[t]-input_knots[t-1];
        int starting_line_index=(t-1)*(order);
        for (int j = 0; j < order; j++)
        {
            add_diff(matrix,starting_line_index+j,order,input_knots,t,j,1);
            matrix(starting_line_index+j,(order+1)*t+j)=-get_fas(j);
        }
    }
    if(!dots.empty()){
        for (int i=0;i < dots.size(); i++) {
            line_index++;
            double value_dimension=value[i];
            int dot=dots[i],difforders=difforder[i];
            if(dot<0||dot>(int)N-1){
                cout<<"no such point:"<<dot<<endl;
                continue;
            }
            if(difforders<0||difforders>order-1){
                cout<<"The order of the derivative is too high:"<<difforders<<endl;
                continue;
            }
            if(dot<N-1) matrix(line_index,dot*(order+1)+difforders)=get_fas(difforders);
            else{ add_diff(matrix,line_index,order,input_knots,dot,difforders,1);}
            target(line_index) = value_dimension;
        }
    }
    if(!dots1.empty()){
        for (int i=0; i < dots1.size(); i++) {
            line_index++;
            int dot1=dots1[i],difforders1=difforder1[i],dot2=dots2[i],difforders2=difforder2[i];
            if(dot1<0||dot1>(int)N-1||dot2<0||dot2>(int)N-1){
                cout<<"no such point: "<<dot1<<" or "<<dot2<<endl;
                continue;
            }
            if(difforders1<0||difforders1>order-1||difforders2<0||difforders2>order-1){
                cout<<"The order of the derivative is too high: "<<difforders1<<" or "<<difforders2<<endl;
                continue;
            }
            if(dot1<N-1) matrix(line_index,dot1*(order+1)+difforders1)=get_fas(difforders1);
            else{ add_diff(matrix,line_index,order,input_knots,dot1,difforders1,1);}
            if(dot2<N-1) matrix(line_index,dot2*(order+1)+difforders2)=-get_fas(difforders2);
            else{ add_diff(matrix,line_index,order,input_knots,dot2,difforders2,-1);}


        }
    }
    if(!exist.empty()){
        for(int i=0;i<(int)exist.size();i++){
            int number=exist[i];
            line_index++;
            if(number<1||number>(int)N-2){
                cout<<"The "<<order<<"th derivative cannot exist at both ends, or out of range "<<endl;
                continue;
            }
            matrix(line_index,number*(order+1)+order)=1;
            matrix(line_index,(number-1)*(order+1)+order)=-1;
            //右n阶导数为以该点起始区间n次系数乘上n!
            //左n阶导数为商议区间n次系数乘上n!
        }

    }
  
    
    process_lines(matrix,target);
    //cout<<matrix<<endl;

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
    if(lu_decomp.rank()<result_size){
        
        cout<<"There isn't unique solution,rank:"<<lu_decomp.rank()<<endl;
        cout<<"It may be due to the choice of conditions that leads to near-singularity."<<endl;
    }
    Eigen::VectorXd solution = matrix.colPivHouseholderQr().solve(target);
    //Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower|Eigen::Upper> cg;
    //cg.compute(matrix);
    //Eigen::VectorXd solution = cg.solve(target);
    

    vector<vector<double>> result;
    for (size_t i=0;i<N-1;i++){
        vector<double> temp;
        for(int j=0;j<=order;j++) temp.push_back(solution[i*(order+1)+j]);
        result.push_back(temp);
    }
    //cout<<target<<endl;
    //cout<<solution<<endl;
    

    return result;


}