#ifndef pputils
#define pputils
#include"../../include/Eigen/Dense"
#include <vector>
void process_lines(Eigen::MatrixXd& A, Eigen::VectorXd& b) {
    int m = A.rows();  
    int n = A.cols(); 

    for (int k = 0; k <= std::min(m, n) - 1; ++k) {
        // 找到当前列最大元素
        int maxRow = k;
        for (int i = k + 1; i < m; ++i) {
            if (fabs(A(i, k)) > fabs(A(maxRow, k))) {
                maxRow = i;
            }
        }
        if (maxRow != k) {
            A.row(k).swap(A.row(maxRow));
            std::swap(b(k), b(maxRow));
        }
        double max_obj=0;
        for (int i = 0; i < n; i++)
        {
            if(fabs(A(k,i))>max_obj) max_obj=fabs(A(k,i));
        }
        if(max_obj>0){
            A.row(k)/=max_obj;
            b[k]/=max_obj;
        }
    }
}


double get_fas(int n){
    if(n<0) throw "invalid";
    if(n==0) return 1;
    else return n*get_fas(n-1);
}




#endif