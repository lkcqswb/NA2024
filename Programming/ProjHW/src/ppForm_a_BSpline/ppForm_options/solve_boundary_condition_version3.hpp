#include<iostream>
#include<vector>
#include<cmath>
#include"../../include/Eigen/Dense"
#include"../../include/Eigen/IterativeLinearSolvers"
#include"ppform_utils.hpp"
using namespace std;
/*
    这里的大致思路为假设最左端结点的1到n-1阶右导数，然后可以利用牛顿插值法递推出所有结点的1到n-1导数。
    结合边界条件的n-1个等式，计算出结果。
    所有结点的1到n-1阶导数均可以表示为最左端结点的1到n-1阶右导数的线性组合。只需要用vector记录系数即可
    构造divided difference
    x f
    t f(t)
    t f(t) f'(t)/1
    ......
*/



Eigen::MatrixXd update_coefficent_new(vector<vector<long long>>C,vector<double>knots,Eigen::MatrixXd co_pre,int index,int order,int size){
    double delta=knots[index]-knots[index-1],delta2=knots[index+1]-knots[index];
    Eigen::VectorXd c0(size);
    c0.setZero();
    if(index!=(int)knots.size()-1){
        c0[index+1]=(double)1/delta2;
        c0[index]=-c0[index+1];
    }
    Eigen::VectorXd c1(size);
    c1.setZero();
    c1[index]=1;

    Eigen::MatrixXd coefficients = co_pre.block(1, 0, co_pre.rows() - 1, co_pre.cols());
    Eigen::MatrixXd matrix1=Eigen::MatrixXd::Zero(order, order);
 
    for (int i = 0; i <order-1; i++)
    {
        for (int j = i; j < order; j++)
        {
            matrix1(i,j)=C[j+1][i+1]*pow(delta,j-i);
        }
    }

    for (int i = 0; i < order-1; i++)
    {
        matrix1(order-1,i)=-pow(delta2,-order)*(pow(delta+delta2,i+1)-pow(delta,i+1));
    }
    matrix1(order-1,order-1)=-pow(delta2,-order)*(pow(delta+delta2,order)-pow(delta,order)-pow(delta2,order));
    Eigen::MatrixXd result=matrix1*coefficients;

    result.row(order-1)+=c0/pow(delta2,order-1);

    Eigen::MatrixXd newMat(result.rows() + 1, result.cols());
    newMat.row(0) = c1;
    newMat.block(1, 0, result.rows(), result.cols()) = result;

    return newMat;
}

vector<vector<long long>> generateMatrix(int n) {
    vector<vector<long long>> C(n + 1, vector<long long>(n + 1, 0));
    for (int i = 0; i <= n; ++i) {
        C[i][0] = 1; // C(n, 0) = 1
        C[i][i] = 1; // C(n, n) = 1
    }
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j < i; ++j) {
            C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
        }
    }
    return C;
}


vector<vector<double>> pp_solve_3(int order,vector<double> input_knots,vector<int> dots1,vector<int> difforder1,vector<int> dots2,vector<int> difforder2,vector<int> dots,vector<int> difforder,vector<double> value,vector<int> exist){
    vector<double> knots=input_knots;
    int N=knots.size();

    knots.push_back(knots[N-1]+1);//用于计算

    int condition_numbers=dots1.size()+dots.size()+exist.size();//共额外条件数
    int result_size=N+order-1;
    if(condition_numbers<result_size){
        cout<<"enadequate information"<<endl;
        throw "enadequate information";
    }
   
    vector<vector<long long>>C= generateMatrix(order);

    vector<Eigen::MatrixXd> coeff;//存储(t1,t2)-(t_{n-1},t_{n})对应点的0-(n)次系数

    //张量运算，每行元素为(f_0,f_1...f_{N-1},m_1/1,m_2/2..m_{n-1}/(n-1)!)对应系数。 m_i为第一个结点i阶导数值，f_i为第i个节点函数值。
    Eigen::MatrixXd t_0=Eigen::MatrixXd::Zero(order+1,result_size);//(t1处的各个系数值)
    t_0(0,0)=1;
    for (int i = 1; i < order; i++) t_0(i,N-1+i)=1;
    t_0(order,1)=1;
    for (int i = 0; i < order; i++) t_0.row(order)=(t_0.row(order)-t_0.row(i))/(knots[1]-knots[0]);
    coeff.push_back(t_0);
   
    for (int i = 1; i <= N-1; i++)//一个个区间更新下去
    {
        coeff.push_back(update_coefficent_new(C,knots,coeff[i-1],i,order,result_size));//这个区间的系数计算下个区间的系数
    }
    


    //构建方程
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(condition_numbers, result_size);
    Eigen::VectorXd target(condition_numbers);

    size_t line_index=0;
    if(!dots.empty()){
        int i=0;
        for (;line_index < dots.size(); line_index++) {
            double value_dimension=value[i];
            int dot=dots[i],difforders=difforder[i];
            if(dot<0||dot>(int)N-1){
                cout<<"no such point:"<<dot<<endl;
                i++;
                continue;
            }
            if(difforders<0||difforders>order-1){
                cout<<"The order of the derivative is too high:"<<difforders<<endl;
                i++;
                continue;
            }
            matrix.row(line_index)=coeff[dot].row(difforders)*get_fas(difforders);
            target(line_index) = value_dimension;
            i++;
        }
    }
    if(!dots1.empty()){
        int i=0;
        for (; line_index < dots1.size()+dots.size(); line_index++) {
            int dot1=dots1[i],difforders1=difforder1[i],dot2=dots2[i],difforders2=difforder2[i];
            if(dot1<0||dot1>(int)N-1||dot2<0||dot2>(int)N-1){
                cout<<"no such point: "<<dot1<<" or "<<dot2<<endl;
                i++;
                continue;
            }
            if(difforders1<0||difforders1>order-1||difforders2<0||difforders2>order-1){
                cout<<"The order of the derivative is too high: "<<difforders1<<" or "<<difforders2<<endl;
                i++;
                continue;
            }
            matrix.row(line_index)+=coeff[dot1].row(difforders1);
            matrix.row(line_index)-=coeff[dot2].row(difforders2);
            target(line_index) =0;
            i++;
        }
    }
    if(!exist.empty()){
        int i=0;
        for(;(int)line_index<condition_numbers;line_index++){
            int number=exist[i];
            if(number<1||number>(int)N-2){
                cout<<"The "<<order<<"th derivative cannot exist at both ends, or out of range "<<endl;
                i++;
                continue;
            }
            matrix.row(line_index)+=coeff[number].row(order);
            matrix.row(line_index)-=coeff[number-1].row(order);
            //右n阶导数为以该点起始区间n次系数乘上n!
            //左n阶导数为商议区间n次系数乘上n!
            target(line_index) = 0;
            i++;
        }

    }
  
    
    process_lines(matrix,target);
    cout<<matrix;

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
    for (size_t i=0;i<coeff.size();i++){
        vector<double> temp;
        for(int j=0;j<=order;j++) temp.push_back(coeff[i].row(j).dot(solution));
        result.push_back(temp);
    }
    
    
    

    return result;


}