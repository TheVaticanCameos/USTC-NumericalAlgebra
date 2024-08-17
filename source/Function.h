#pragma once
#include <iostream>
#include <vector>
using namespace std;

// =======================第一章中声明的函数===================

void forward_subs(vector<vector<double>> L, vector<double>& b);//前代法

void forward_subs1(vector<vector<double>> L, vector<double>& b);//对角元为1的前代法

void back_subs(vector<vector<double>> U, vector<double>& b);//回代法

void back_subs1(vector<vector<double>> U, vector<double>& y);//对角元为1的回代法

void gauss_elim(vector<vector<double>>& A);//Gauss消去法

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//全主元Gauss消去法

void gauss_elim_col_pivoting( vector<vector<double>>& A, vector<int>& u);//列主元Gauss消去法

void vector_pb(vector<int>&u,vector<double>&b);//计算向量P*b【可选】

void vector_qb(vector<int>& v, vector<double>& b);//计算向量Q*b【可选】

void cholesky_decomp(vector<vector<double>>& A);//对称正定阵标准Cholesky分解

void modified_cholesky_decomp(vector<vector<double>>& A);//改进的平方根法

void matrix_DLT(vector<vector<double>>& A);//计算矩阵D*L^T【可选】



// =======================第二章中声明的函数===================

void mat_mul_vec(vector<vector<double> > A, vector<double> x, vector<double>& b);	// 矩阵和向量的乘法

void sol_eqs_col_LU(vector<vector<double> > A, vector<double> b, vector<double>& sol);	// 列主元Gauss消去法求解线性方程组

double vec_1_norm(vector<double> w);	// 向量的1范数

double vec_infty_norm(vector<double> z, int& j);	// 向量的无穷范数

double vec_infty_norm(vector<double> z);	// 计算向量z的无穷范数（重载）

double vec_inner_prod(vector<double> x, vector<double> y);	// 两个向量的内积

void transpose(vector<vector<double> > A, vector<vector<double> >& At);	// 矩阵的转置

void vec_sign(vector<double> w, vector<double>& v);	// 向量的符号函数

double inv_mat_infty_norm(vector<vector<double> > A);	// 逆矩阵的无穷范数估计（优化法）

double mat_infty_norm(vector<vector<double> > A);	// 矩阵的无穷范数

double kappa_infty(vector<vector<double> > A);	// 无穷范数意义下的条件数

void print_mat(vector<vector<double> > A);	// 打印矩阵

void print_vec(vector<double> b);	// 打印向量

double precision(vector<vector<double> > A, vector<double> b, vector<double> x);	// 计算解的精度估计



// =======================第三章中声明的函数===================

double householder(vector<double> x, vector<double>& v);	// 计算Householder变换，对应于算法3.2.1

vector<double> mat_slice_vec(vector<vector<double> > A, int i, int j, int k);	// 给定矩阵A，返回向量A(i:j, k)

vector<vector<double> > mat_slice_mat(vector<vector<double> > A, int i1, int i2, int j1, int j2);	// 给定矩阵A，返回子矩阵A(i1:i2, j1:j2)

vector<vector<double> > mat_sub(vector<vector<double> > A, vector<vector<double> > B);	// 计算矩阵减法A-B，返回相减得到的矩阵

vector<vector<double> > mat_I_sub_beta_vvT(double beta, vector<double> v);	// 计算I-beta*v*v^T，返回计算结果(矩阵)

vector<vector<double> > mat_mul_mat(vector<vector<double> > A, vector<vector<double> > B);	// 计算矩阵乘法A*B，返回计算结果

double vec_2_norm(vector<double> x);	// 计算向量x的2范数

vector<double> QR_decomp_householder(vector<vector<double> >& A);	// 用Householder变换实现矩阵A的QR分解，对应于算法3.3.1

void QR_solver(vector<vector<double> > A, vector<double>& b);	// 利用QR分解，求解线性方程组Ax=b

double QR_LS(vector<vector<double> > A, vector<double> b, vector<double>& x);	// 利用QR分解，求解最小二乘问题 min ||Ax-b||_2



// =======================第四章中声明的函数===================

vector<double> Jacobi_iteration(vector<double> x, vector<vector<double> > B, vector<double> g);	// Jacobi迭代法，计算迭代1步得到的向量

void Jacobi_preparation(vector<vector<double> >& A, vector<double>& b);	// 计算Jacobi迭代的系数矩阵B、常向量g

vector<double> Gauss_Seidel_iteration(vector<double> x, vector<vector<double> > B, vector<double> g);	// Gauss-Seidel迭代法，计算迭代1步得到的向量

vector<double> SOR_iteration(vector<double> x, vector<vector<double> > A, vector<double> b, double omega);	// SOR迭代法，计算1步迭代得到的向量

vector<double> vec_sub(vector<double> a, vector<double> b);	// 计算两个向量a, b的差，即a-b，返回计算结果

double mat_flatten_2_norm(vector<vector<double> > A);	// 计算矩阵A展平为向量的2范数

vector<vector<double> > pde2D_Jacobi(int n, double eps);	// Jacobi迭代法，五点差分格式

double find_min(vector<vector<double> > A, int& ii, int& jj);	// 找到解中的最小分量并返回，取到最小值的下标为(ii, jj)

vector<vector<double> > pde2D_GS(int n, double eps);	// Gauss-Seidel迭代法，五点差分格式

vector<vector<double> > pde2D_SOR(int n, double eps, double omega);	// SOR迭代法，五点差分格式



// =======================第五章中声明的函数===================

void CGD(vector<vector<double> > A, vector<double>& x, vector<double> b, double eps, int kmax);	// 共轭梯度法求解线性方程组Ax=b



// =======================第六章中声明的函数===================

double powerIterationMaxEigenvalue(vector<vector<double> > A, vector<double> u);	// 幂法求矩阵A的模最大的特征值，返回计算结果

double maxAbsRoot(vector<double> coef);	// 使用幂法，计算并返回多项式方程的模最大的根

vector<vector<double> > eye(int n);	// 生成一个n阶单位阵并返回

void mat_put_mat(vector<vector<double> >& A, vector<vector<double> > subMat, int i, int j);	// 将矩阵subMat放入矩阵A中，其中subMat[0][0]放在A[i][j]的位置

void upperHessenberg(vector<vector<double> >& A);	// 矩阵A的上Hessenberg化

void doubleShiftQR(vector<vector<double> >& H, vector<vector<double> >& P);	// Francis双重步位移的QR迭代算法，对应于教材算法6.4.2

vector<vector<double> > transpose(vector<vector<double> > A);	// 返回矩阵A的转置

bool isApproximable(vector<vector<double> > H);	// 判断Hessenberg矩阵H是否可约

bool isProposedUpTriMat(vector<vector<double> > H);	// 判断Hessenberg矩阵H是否是拟上三角阵

void find_m_l(vector<vector<double> > H, int& m, int& l);	// 算法6.4.3的辅助函数，完成(3)(ii)的功能，即找到满足要求的最大的非负整数m和最小的非负整数l

void implicitQR(vector<vector<double> >& A, double u);	// 计算实矩阵A的实Schur分解：隐式QR算法，对应于教材算法6.4.3

void implicitQREigenvalue(vector<vector<double> > A, vector<double>& Re, vector<double>& Im);	// 使用隐式QR迭代，计算矩阵A的全部特征值



// =======================第七章中声明的函数===================

double nonDiagNorm(vector<vector<double> > A);	// 计算方阵A的非对角“范数”，其定义见教材P211

bool isPassed(vector<vector<double> > A, double delta, int& p, int& q);	// 判断过关Jacobi方法中，当前矩阵是否“过关”

void rotateJacobi(vector<vector<double> > A, int p, int q, double& c, double& s);	// 计算J(p,q,\theta)中的\cos\theta和\sin\theta，保存在c, s中

vector<vector<double> > leftJacobiMul(vector<vector<double> > A, int p, int q, double c, double s);	// 用Jacobi矩阵J(p,q,\theta)左乘矩阵A，返回计算结果

vector<vector<double> > rightJacobiMul(vector<vector<double> > A, int p, int q, double c, double s);	// 用Jacobi矩阵J(p,q,\theta)右乘矩阵A，返回计算结果

void passingJacobiMethod(vector<vector<double> > A, vector<double>& EigenValues, vector<vector<double> >& EigenVectors, double sigma);  // 过关Jacobi方法求对称矩阵的全部特征值和特征向量

void print_mat(vector<vector<double> > A, int digits);	// 打印矩阵A，其中digits为小数点后位数（重载）

void print_vec(vector<double> v, int digits);	// 打印向量v，其中digits为小数点后位数（重载）

int changeSignNum(vector<vector<double> > A, double mu);	// 计算 对称三对角矩阵 A在mu处的变号数s_n(mu)

double bisectionMethod(vector<vector<double> > A, int m, double eps);	// 使用二分法，求矩阵A的第m小的特征值

vector<double> inversePowerMethod(vector<vector<double> > A, double lambda, double eps);	// 使用反幂法求解矩阵A的特征值lambda对应的特征向量

void helper(vector<double> eigen);	// 打表格小助手



// =======================第八章中声明的函数===================

vector<vector<double> > xyT(vector<double> x, vector<double> y);	// 计算矩阵 xy^T，其中x, y为行向量，大小分别为m, n

vector<double> numMul(double lambda, vector<double> v);	// 计算向量数乘，返回计算结果

vector<vector<double> > numMul(double lambda, vector<vector<double> > A);	// 计算矩阵数乘，返回计算结果

vector<vector<double> > biDiagonalization(vector<vector<double> > A, vector<vector<double> >& U, vector<vector<double> >& V);	// 对矩阵 A 用 Householder 变换做二对角化

void givens(double& c, double& s, double a, double b);	// 找Givens变换中的\cos\theta和\sin\theta，分别保存在c, s中

vector<vector<double> > WilkinsonSVD(vector<vector<double> > B, vector<vector<double> >& P, vector<vector<double> >& Q);	// 带Wilkinson位移的SVD迭代

void find_p_q(vector<vector<double> > B, int& p, int& q);	// 算法7.6.3中的辅助程序，找满足要求的p和q

void SVD(vector<vector<double> > A, vector<vector<double> >& U, vector<vector<double> >& Sigma, vector<vector<double> >& V, double eps);	// 矩阵A的SVD分解，其中A的大小为m*n，m\geq n

double maxAbs(vector<vector<double> > A);	// 返回矩阵A中模最大的元素