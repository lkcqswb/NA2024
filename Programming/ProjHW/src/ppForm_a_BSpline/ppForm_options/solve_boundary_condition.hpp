#include<iostream>
#include<vector>
#include<cmath>
#include"../../include/Eigen/Dense"
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

double vec_get_solution(vector<double> x,Eigen::VectorXd solution,int order){
    if((int)(x.size())!=order||(int)(solution.size())!=order-1){
        cout<<"invalid shape"<<endl;
        throw "invalid shape";
    }
    return x[0] + (Eigen::Map<const Eigen::VectorXd>(x.data() + 1, order - 1).dot(solution));
}


vector<vector<double>> update_deri(vector<double>knots,vector<vector<double>> co_previous,int index,int order){
    vector<vector<double>> t_now;
    vector<vector<double>> co_pre=co_previous;
    for (int j = 1; j < order; j++)//从1阶导数到n-1阶导数。
        {
            vector<double> t_i(order,0);
            for (int k = j; k <=order; k++)
            {
                co_pre[k-1]=num_product(co_pre[k-1],(k+1-j));//对每一项进行求导。
            }
            for (int k = j-1; k < order; k++) t_i=vec_add(t_i,num_product(co_pre[k],pow(knots[index]-knots[index-1],k+1-j)));
            t_now.push_back(t_i);
    }
    return t_now;
}

vector<vector<double>> update_coefficent(vector<double> factorial_index,vector<double>knots,vector<vector<double>> t_now,vector<double> f_values,int index,int order){
    vector<vector<double>> co_now;
    for (int j = 1; j < order; j++)//从第一项到n-1项系数。
    {
        co_now.push_back(num_product(t_now[j-1],(double)(1)/factorial_index[j]));
    }
    vector<double> last_coef(order,0);//第n项系数需要单独计算
    last_coef[0]=(f_values[index+1]-f_values[index])/(knots[index+1]-knots[index]);
    for (int i = 1; i < order; i++) last_coef=num_product(vec_add(last_coef,vec_minus(co_now[i-1])),1/(knots[index+1]-knots[index]));
    co_now.push_back(last_coef);
    return co_now;
}




vector<vector<double>> pp_solve(int order,vector<double> knots,vector<double> f_values,vector<int> dots1,vector<int> difforder1,vector<int> dots2,vector<int> difforder2,vector<int> dots,vector<int> difforder,vector<double> value,vector<int> exist){

    int condition_numbers=dots1.size()+dots.size()+exist.size();//共额外条件数

    if(condition_numbers<order-1){
        cout<<"enadequate information"<<endl;
        throw "enadequate information";
    }
    if(order==1){
        vector<vector<double>> result;
        for (size_t i=0;i<knots.size()-1;i++){
            result.push_back({f_values[i],(f_values[i+1]-f_values[i])/(knots[i+1]-knots[i])});
        }
        return result;
    }

    vector<double> factorial_index(1,1);
    for (int i = 1; i < order; i++) factorial_index.push_back(i*factorial_index[i-1]);//存储1到n-1的阶乘，免于后续反复计算。



    vector<vector<vector<double>>> deriv_co;//存储t1-tn对应点的1-(n-1)阶导数
    vector<vector<vector<double>>> coeff;//存储(t1,t2)-(t_{n-1},t_{n})对应点的1-(n)次系数


    vector<vector<double>> t_0;;//(t1处的1阶导数到n-1阶导数)
    for (int i = 1; i < order; i++)
    {
        vector<double> deri(order,0);//i阶导数
        deri[i]=1;
        t_0.push_back(deri);
    }
    deriv_co.push_back(t_0);
    coeff.push_back(update_coefficent(factorial_index,knots,t_0,f_values,0,order));//计算第一个区间对应的系数

    for (size_t i = 1; i < knots.size()-1; i++)//一个个区间更新下去
    {
        deriv_co.push_back(update_deri(knots,coeff[i-1],(int)i,order));//上一个区间的系数，用于计算这个区间右侧点的导数
        coeff.push_back(update_coefficent(factorial_index,knots,deriv_co[i],f_values,i,order));//这个区间的导数计算下个区间的系数
    }
    deriv_co.push_back(update_deri(knots,coeff[coeff.size()-1],(int)(knots.size()-1),order));//计算最右端点的导数


    //构建方程
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(condition_numbers, order - 1);
    Eigen::VectorXd target(condition_numbers);

    size_t line_index=0;
    if(!dots.empty()){
        int i=0;
        for (;line_index < dots.size(); line_index++) {
            double value_dimension=value[i];
            int dot=dots[i],difforders=difforder[i];
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
            vector<double> derivation = deriv_co[dot][difforders - 1];
            for (int j = 1; j < order; j++) {
                matrix(line_index, j - 1) = derivation[j];
            }
            target(line_index) = value_dimension - derivation[0];
            i++;
        }
    }
    if(!dots1.empty()){
        int i=0;
        for (; line_index < dots1.size()+dots.size(); line_index++) {
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
            vector<double> derivation1 = deriv_co[dot1][difforders1 - 1];
            vector<double> derivation2 = deriv_co[dot2][difforders2 - 1];
            for (int j = 1; j < (int)derivation1.size(); j++) matrix(line_index, j-1) += derivation1[j];
            for (int j = 1; j < (int)derivation2.size(); j++) matrix(line_index, j-1) -= derivation2[j];
            target(line_index) =derivation2[0]-derivation1[0] ;
            i++;
        }
    }
    if(!exist.empty()){
        int i=0;
        for(;(int)line_index<condition_numbers;line_index++){
            int number=exist[i];
            if(number<1||number>(int)knots.size()-2){
                cout<<"The "<<order<<"th derivative cannot exist at both ends, or out of range "<<endl;
                continue;
            }
            vector<double> derivation1 =coeff[number][order-1];//右n阶导数为以该点起始区间n次系数乘上n!
            vector<double> derivation2 =coeff[number-1][order-1];//左n阶导数为商议区间n次系数乘上n!
            for (int j = 1; j < (int)derivation1.size(); j++) matrix(line_index, j-1) += derivation1[j];
            for (int j = 1; j < (int)derivation2.size(); j++) matrix(line_index, j-1) -= derivation2[j];
            target(line_index) = derivation2[0]-derivation1[0];
            i++;
        }

    }









    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
    if(lu_decomp.rank()<order-1){
        cout<<"There isn't unique solution"<<endl;
        throw "There isn't unique solution";
    }
    Eigen::VectorXd solution = matrix.colPivHouseholderQr().solve(target);
    
    

    vector<vector<double>> result;
    for (size_t i=0;i<coeff.size();i++){
        vector<double> temp={f_values[i]};
        for(int j=0;j<order;j++) temp.push_back(vec_get_solution(coeff[i][j],solution,order));
        result.push_back(temp);
    }

    

    return result;


}