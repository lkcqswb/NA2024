#include<iostream>
#include<vector>
#include<cmath>
#include"../../include/Eigen/Dense"
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


vector<double> num_product(vector<double> x,double b){
    vector<double> result(x.size());
    transform(x.begin(), x.end(), result.begin(), [b](double val) { return val * b; });
    return result;
}

vector<double> vec_add(vector<double> x,vector<double> y){
    if(x.size()!=y.size()){
        cout<<"illegal add"<<endl;
        throw "illegal add";
    }
    vector<double> result(x.size());
    std::transform(x.begin(), x.end(), y.begin(), result.begin(), plus<double>());
    return result;
}

vector<double> vec_minus(vector<double> x){
    vector<double> result(x.size());
    transform(x.begin(), x.end(), result.begin(), [](double val) { return -val; });
    return result;
}

double vec_get_solution(const std::vector<double>& x, const Eigen::VectorXd& solution) {
    if ((int)x.size() != (int)solution.size() + 1) {
        throw std::invalid_argument("Vectors must have compatible sizes.");
    }

    double firstElement = x[0];
    Eigen::Map<const Eigen::VectorXd> subX(x.data() + 1, x.size() - 1);
    return firstElement + (subX.dot(solution));
}


vector<vector<double>> update_deri(vector<double>knots,vector<vector<double>> co_previous,int index,int order,vector<double> f_value){
    vector<vector<double>> t_now;
    vector<vector<double>> co_pre=co_previous;
    //0阶导数
    vector<double> t_0(order+1,0);
    if(index<(int)knots.size()-1) t_0[0]=f_value[index];
    else t_0[1]=1;
    t_now.push_back(t_0);
    for (int j = 1; j < order; j++)//从1阶导数到n-1阶导数。
        {
            vector<double> t_i(order+1,0);
            for (int k = j; k <=order; k++)
            {
                co_pre[k]=num_product(co_pre[k],(k+1-j));//对每一项进行求导。
            }
            for (int k = j; k <=order; k++) t_i=vec_add(t_i,num_product(co_pre[k],pow(knots[index]-knots[index-1],k-j)));
            t_now.push_back(t_i);
    }
    return t_now;
}

vector<vector<double>> update_coefficent(vector<double> factorial_index,vector<double>knots,vector<vector<double>> t_now,int index,int order,vector<double> f_value){
    vector<vector<double>> co_now;
    for (int j = 0; j < order; j++)//从第0次项到n-1次项系数。
    {
        co_now.push_back(num_product(t_now[j],(double)(1)/factorial_index[j]));
    }
    vector<double> last_coef(order+1,0);//第n项系数需要单独计算
    if(index+1<(int)knots.size()-1) last_coef[0]=f_value[index+1];
    else last_coef[1]=1;

    for (int i = 0; i < order; i++) last_coef=num_product(vec_add(last_coef,vec_minus(co_now[i])),1/(knots[index+1]-knots[index]));
    co_now.push_back(last_coef);
    return co_now;
}




vector<vector<double>> pp_solve_2(int order,vector<double> knots,vector<int> dots1,vector<int> difforder1,vector<int> dots2,vector<int> difforder2,vector<int> dots,vector<int> difforder,vector<double> value,vector<int> exist){
    vector<double> f_values(knots.size()-1,0);//提取0-n-2点的函数值。
    vector<double> new_dots,new_difforders,new_value;
    if(!dots.empty()){
        for (size_t i=0;i<dots.size();i++) {
            double value_dimension=value[i];
            int dot=dots[i],difforders=difforder[i];
            if(difforders==0&&dot>=0&&dot<=(int)knots.size()-2){
                f_values[dot]=value_dimension;
            }else{
                new_dots.push_back(dots[i]);
                new_difforders.push_back(difforder[i]);
                new_value.push_back(value[i]);
            }
        }
    }
    int condition_numbers=dots1.size()+new_dots.size()+exist.size();//共额外条件数

    if(condition_numbers<order-1){
        cout<<"enadequate information"<<endl;
        throw "enadequate information";
    }

    vector<double> factorial_index(1,1);
    for (int i = 1; i < order; i++) factorial_index.push_back(i*factorial_index[i-1]);//存储1到n-1的阶乘，免于后续反复计算。

    vector<vector<vector<double>>> deriv_co;//存储t1-tn对应点的0-(n-1)阶导数
    vector<vector<vector<double>>> coeff;//存储(t1,t2)-(t_{n-1},t_{n})对应点的0-(n)次系数

    //张量运算，设置每个元素为(常数项,f_{N-1},m_1,m_2..m_{n-1})对应系数。 m_i为第一个结点i阶导数值，f_i为第i个节点函数值。
    vector<vector<double>> t_0;;//(t1处的0阶导数到n-1阶导数)
    for (int i = 0; i < order; i++)
    {
        vector<double> deri(order+1,0);//i阶导数
        if(i==0) deri[0]=f_values[0];
        else deri[1+i]=1;
        t_0.push_back(deri);
    }
    
    deriv_co.push_back(t_0);
    coeff.push_back(update_coefficent(factorial_index,knots,t_0,0,order,f_values));//计算第一个区间对应的系数

    for (size_t i = 1; i < knots.size()-1; i++)//一个个区间更新下去
    {
        deriv_co.push_back(update_deri(knots,coeff[i-1],(int)i,order,f_values));//上一个区间的系数，用于计算这个区间右侧点的导数
        coeff.push_back(update_coefficent(factorial_index,knots,deriv_co[i],i,order,f_values));//这个区间的导数计算下个区间的系数
    }
    deriv_co.push_back(update_deri(knots,coeff[coeff.size()-1],(int)(knots.size()-1),order,f_values));//计算最右端点的导数


    //构建方程
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(condition_numbers, order);
    Eigen::VectorXd target(condition_numbers);

    size_t line_index=0;
    if(!new_dots.empty()){
        int i=0;
        for (;line_index < new_dots.size(); line_index++) {
            double value_dimension=new_value[i];
            int dot=new_dots[i],difforders=new_difforders[i];
            if(dot<0||dot>(int)knots.size()-1){
                cout<<"no such point:"<<dot<<endl;
                i++;
                continue;
            }
            if(difforders<0||difforders>order-1){
                cout<<"The order of the derivative is too high:"<<difforders<<endl;
                i++;
                continue;
            }
            vector<double> derivation = deriv_co[dot][difforders];
            for (int j = 1; j < (int)derivation.size(); j++) {
                matrix(line_index, j-1) = derivation[j];
            }
            target(line_index) = value_dimension-derivation[0];
            i++;
        }
    }
    if(!dots1.empty()){
        int i=0;
        for (; line_index < dots1.size()+new_dots.size(); line_index++) {
            int dot1=dots1[i],difforders1=difforder1[i],dot2=dots2[i],difforders2=difforder2[i];
            if(dot1<0||dot1>(int)knots.size()-1||dot2<0||dot2>(int)knots.size()-1){
                cout<<"no such point: "<<dot1<<" or "<<dot2<<endl;
                i++;
                continue;
            }
            if(difforders1<0||difforders1>order-1||difforders2<0||difforders2>order-1){
                cout<<"The order of the derivative is too high: "<<difforders1<<" or "<<difforders2<<endl;
                i++;
                continue;
            }
            vector<double> derivation1 = deriv_co[dot1][difforders1];
            vector<double> derivation2 = deriv_co[dot2][difforders2];
            for (int j = 1; j < (int)derivation1.size(); j++) matrix(line_index, j-1) += derivation1[j];
            for (int j = 1; j < (int)derivation2.size(); j++) matrix(line_index, j-1) -= derivation2[j];
            target(line_index) =derivation2[0]-derivation1[0];
            i++;
        }
    }
    if(!exist.empty()){
        int i=0;
        for(;(int)line_index<condition_numbers;line_index++){
            int number=exist[i];
            if(number<1||number>(int)knots.size()-2){
                cout<<"The "<<order<<"th derivative cannot exist at both ends, or out of range "<<endl;
                i++;
                continue;
            }
            vector<double> derivation1 =coeff[number][order];//右n阶导数为以该点起始区间n次系数乘上n!
            vector<double> derivation2 =coeff[number-1][order];//左n阶导数为商议区间n次系数乘上n!
            for (int j = 1; j < (int)derivation1.size(); j++) matrix(line_index, j-1) += derivation1[j];
            for (int j = 1; j < (int)derivation2.size(); j++) matrix(line_index, j-1) -= derivation2[j];
            target(line_index) = derivation2[0]-derivation1[0];
            i++;
        }

    }

    
    

    process_lines(matrix,target);
    
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
    if((int)lu_decomp.rank()<order){
        
        cout<<"There isn't unique solution,rank:"<<lu_decomp.rank()<<endl;
        cout<<"It may be due to the choice of conditions that leads to near-singularity."<<endl;
    }
    Eigen::VectorXd solution = matrix.colPivHouseholderQr().solve(target);

    

    vector<vector<double>> result;
    for (size_t i=0;i<coeff.size();i++){
        vector<double> temp;
        for(int j=0;j<=order;j++) temp.push_back(vec_get_solution(coeff[i][j],solution));
        result.push_back(temp);
    }

    return result;

}