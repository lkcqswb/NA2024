#include<cmath>
#include<iostream>
#include"../include/Eigen/Dense"
using namespace std;

void gaussianEliminationWithPartialPivoting(Eigen::MatrixXd& A, Eigen::VectorXd& b) {
    int m = A.rows();  
    int n = A.cols();  

    for (int k = 0; k < std::min(m, n) - 1; ++k) {
        // 找到当前列的主元
        int maxRow = k;
        for (int i = k + 1; i < m; ++i) {
            if (fabs(A(i, k)) > fabs(A(maxRow, k))) {
                maxRow = i;
            }
        }

        // 行交换
        if (maxRow != k) {
            A.row(k).swap(A.row(maxRow));
            std::swap(b(k), b(maxRow));
        }

        // 消元
        for (int i = k + 1; i < m; ++i) {
            if (A(k, k) != 0) {  // 防止除以零
                double factor = A(i, k) / A(k, k);
                A.row(i) -= factor * A.row(k);
                b(i) -= factor * b(k);
            }
        }
    }
}

int main(){//z^=x^2+y^2
    double t;
    int n;
    cin>>n;
    Eigen::VectorXd vec=Eigen::VectorXd::Zero(n);
    Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(n, n);
    for (size_t i = 0; i < n-2; i++)
    {
        matrix(i,i)=1;
        matrix(i,i+1)=4;
        matrix(i,i+2)=1;
    }
    matrix(n-2,0)=-1;
    matrix(n-2,1)=4;
    matrix(n-2,2)=-6;
    matrix(n-2,3)=4;
    matrix(n-2,4)=-1;
    matrix(n-1,n-5)=-1;
    matrix(n-1,n-4)=4;
    matrix(n-1,n-3)=-6;
    matrix(n-1,n-2)=4;
    matrix(n-1,n-1)=-1;

    //gaussianEliminationWithPartialPivoting(matrix,vec);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double conditionNumber = svd.singularValues()(0) / svd.singularValues().tail<1>()(0);
    std::cout << "Condition number: " << conditionNumber << std::endl;  
    int rank = (svd.singularValues().array() > 1e-10).count();
    cout<<rank;
    double minSingularValue = svd.singularValues().tail<1>()(0);
    std::cout << "Minimum singular value: " << minSingularValue << std::endl;
    
}