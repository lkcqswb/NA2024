#include<iostream>
#include<vector>
#include<string>
#include"../../include/Eigen/Dense"
using namespace std;
double get_ti(int index,vector<double> knots){
    if(index<0) return knots[0]+index;
    if(index<=(int)knots.size()-1) return knots[index];
    else return knots[knots.size()-1]-knots.size()+1+index;
}

int max(int x,int y){
    if(x>y) return x;
    return y;
}
//构造在t[offset]到t[-1]每个点对应
//B^0_0(t_index) B^0_1(t_index) B^0_2(t_index) ... B^0_{index}(t_index) 
//...
//B^order_0(t_index) B^order_1(t_index) B^order_2(t_index) ... B^order_{index}(t_index) 
vector<double> update_B(vector<double> B,vector<double> knots,int current_order,double x){
    for (int k = max(0,B.size()-2-current_order); k <= (int)B.size()-2; k++)
    {
        B[k]=B[k]*(x-get_ti(k-1,knots))/(get_ti(current_order+k,knots)-get_ti(k-1,knots))+B[k+1]*(get_ti(k+current_order+1,knots)-x)/(get_ti(current_order+k+1,knots)-get_ti(k,knots));
    }
    B[B.size()-1]*=(x-get_ti(B.size()-2,knots))/(get_ti(current_order+B.size()-1,knots)-get_ti(B.size()-2,knots));
    return B;
}


vector<vector<double>> construct_value_table(vector<double> knots,int order,int index,double x){
    vector<vector<double>> table={};
    vector<double> B(index+1,0);
    B[index]=1;
    table.push_back(B);
    for (int i = 0; i < order; i++)
    {
        table.push_back(update_B(table[table.size()-1],knots,i,x));
    }
    return table;
}

//构造在t[offset]到t[-1]每个点对应的各个B样条基的1到order-1阶导。
//B^order_0^(n)(t_index) B^order_1(t_index) B^order_2(t_index) ... B^order_{index}(t_index)
//可以用矩阵运算简化,k阶导表示为
// 
//B^n_i^(k)(x)=(B^{n-k}_0^(n)(x),B^{n-k}_1^(n)(x),B^{n-k}_2^(n)(x)...)的线性组合,A最初为对角阵
//再左乘上对应行的value_table就可以得到导数值
vector<vector<double>> construct_derivatives_table(vector<vector<double>> value_table,vector<double> knots){
    vector<vector<double>> difftable;
    int order=value_table.size()-1;
    int index=value_table[0].size()-1;
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(index+1, index+1);
    Eigen::MatrixXd coefficients=Eigen::MatrixXd::Identity(index+1, index+1);
    for (int i = 0; i < order; i++)
    {
        Eigen::MatrixXd values(1,index+1);
        values.row(0)=Eigen::Map<const Eigen::VectorXd>(value_table[order-i-1].data(),index+1);
        for (int j = 0; j < index+1; j++)
        {
            matrix(j,j)=(double)(order-i)/(get_ti(j+order-i-1,knots)-get_ti(j-1,knots));
            if(j<index) matrix(j+1,j)=-(double)(order-i)/(get_ti(j+order-i,knots)-get_ti(j,knots));
        }
        coefficients=matrix*coefficients;
        //cout << "MatrixXd:\n" << coefficients << std::endl;
        Eigen::MatrixXd deriv=values*coefficients;
        vector<double> result(deriv.data(), deriv.data() + deriv.rows() * deriv.cols());
        difftable.push_back(result);
        //cout<<"value"<<values<<endl;
        //cout<<"deriv"<<deriv<<endl;
    }
    //cout << "MatrixXd:\n" << coefficients << std::endl;
    
    return difftable;
}



vector<double> B_solve(int order,vector<double> knots,vector<double> f_values,vector<int> dots1,vector<int> difforder1,vector<int> dots2,vector<int> difforder2,vector<int> dots,vector<int> difforder,vector<double> value,vector<int> exist ){
    vector<vector<vector<double>>> value_table;
    int condition_numbers=dots1.size()+dots.size()+exist.size();//共额外条件数
    

    if(condition_numbers<order-1){
        cout<<"enadequate information"<<endl;
        throw "enadequate information";
    }

    int offset=order-1;//索引order-1开始有对应的函数值
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(condition_numbers+knots.size()-order+1, knots.size());
    Eigen::VectorXd target(condition_numbers+knots.size()-order+1);
    
    
    for (int i = offset; i < (int)knots.size(); i++)
    {
        value_table.push_back(construct_value_table(knots,order,i,knots[i]));

        for(int j=0;j<i+1;j++) matrix(i-offset,j)=value_table[i-offset][order][j];
        target[i-offset]=f_values[i-offset];
    }
    size_t line_index = knots.size()-offset;
    
    if(!dots.empty()){
        int i=0;
        for (;line_index < dots.size()+knots.size()-order+1; line_index++) {
            double value_dimension=value[i];
            int dot=dots[i],difforders=difforder[i];
            if(dot<0||dot>(int)knots.size()-order){
                cout<<"no such point:"<<dot<<endl;
                i++;
                continue;
            }
            if(difforders<0||difforders>order-1){
                cout<<"The order of the derivative is too high:"<<difforders<<endl;
                i++;
                continue;
            }
            vector<vector<double>> derivation =construct_derivatives_table(value_table[dot],knots);
            for (int j = 0; j < (int)derivation[0].size(); j++) {
                matrix(line_index, j) = derivation[difforders-1][j];
            }
            //cout<<"ma:"<<matrix.row(line_index)<<endl;
            target(line_index) = value_dimension;
            i++;
        }
    }
    if(!dots1.empty()){
        int i=0;
        for (; line_index < dots1.size()+dots.size()+knots.size()-order+1; line_index++) {
            int dot1=dots1[i],difforders1=difforder1[i],dot2=dots2[i],difforders2=difforder2[i];
            if(dot1<0||dot1>(int)knots.size()-order||dot2<0||dot2>(int)knots.size()-order){
                cout<<"no such point: "<<dot1<<" or "<<dot2<<endl;
                i++;
                continue;
            }
            if(difforders1<0||difforders1>order-1||difforders2<0||difforders2>order-1){
                cout<<"The order of the derivative is too high: "<<difforders1<<" or "<<difforders2<<endl;
                i++;
                continue;
            }
            vector<vector<double>> derivation1 =construct_derivatives_table(value_table[dot1],knots);
            vector<vector<double>> derivation2 =construct_derivatives_table(value_table[dot2],knots);
            for (int j = 0; j < (int)derivation1[0].size(); j++) matrix(line_index, j) += derivation1[difforders1-1][j];
            for (int j = 0; j < (int)derivation2[0].size(); j++) matrix(line_index, j) -= derivation2[difforders2-1][j];
            target(line_index) = 0;
            i++;
        }
    }
    if(!exist.empty()){
        int i=0;
        for(;line_index<condition_numbers+knots.size()-order+1;line_index++){
            int number=exist[i];
            if(number<1||number>(int)knots.size()-order-1){
                cout<<"The "<<order<<"th derivative cannot exist at both ends, or out of range "<<endl;
                continue;
            }
            vector<vector<double>> derivation1 =construct_derivatives_table(value_table[number],knots);
            vector<vector<double>> derivation2 =construct_derivatives_table(value_table[number+1],knots);
            for (int j = 0; j < (int)derivation1[0].size(); j++) matrix(line_index, j) += derivation1[order-1][j];
            for (int j = 0; j < (int)derivation2[0].size(); j++) matrix(line_index, j) -= derivation2[order-1][j];
            target(line_index) = 0;
            i++;
        }

    }



    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
    if(lu_decomp.rank()<(int)knots.size()){
        cout<<"There isn't unique solution"<<endl;
        throw "There isn't unique solution";
    }
    Eigen::VectorXd solution = matrix.colPivHouseholderQr().solve(target);
    vector<double> result(solution.data(), solution.data() + solution.size());
    //cout<<"sol:"<<solution<<endl;
    return result;
}