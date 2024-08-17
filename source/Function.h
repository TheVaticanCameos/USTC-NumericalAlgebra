#pragma once
#include <iostream>
#include <vector>
using namespace std;

// =======================��һ���������ĺ���===================

void forward_subs(vector<vector<double>> L, vector<double>& b);//ǰ����

void forward_subs1(vector<vector<double>> L, vector<double>& b);//�Խ�ԪΪ1��ǰ����

void back_subs(vector<vector<double>> U, vector<double>& b);//�ش���

void back_subs1(vector<vector<double>> U, vector<double>& y);//�Խ�ԪΪ1�Ļش���

void gauss_elim(vector<vector<double>>& A);//Gauss��ȥ��

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//ȫ��ԪGauss��ȥ��

void gauss_elim_col_pivoting( vector<vector<double>>& A, vector<int>& u);//����ԪGauss��ȥ��

void vector_pb(vector<int>&u,vector<double>&b);//��������P*b����ѡ��

void vector_qb(vector<int>& v, vector<double>& b);//��������Q*b����ѡ��

void cholesky_decomp(vector<vector<double>>& A);//�Գ��������׼Cholesky�ֽ�

void modified_cholesky_decomp(vector<vector<double>>& A);//�Ľ���ƽ������

void matrix_DLT(vector<vector<double>>& A);//�������D*L^T����ѡ��



// =======================�ڶ����������ĺ���===================

void mat_mul_vec(vector<vector<double> > A, vector<double> x, vector<double>& b);	// ����������ĳ˷�

void sol_eqs_col_LU(vector<vector<double> > A, vector<double> b, vector<double>& sol);	// ����ԪGauss��ȥ��������Է�����

double vec_1_norm(vector<double> w);	// ������1����

double vec_infty_norm(vector<double> z, int& j);	// �����������

double vec_infty_norm(vector<double> z);	// ��������z������������أ�

double vec_inner_prod(vector<double> x, vector<double> y);	// �����������ڻ�

void transpose(vector<vector<double> > A, vector<vector<double> >& At);	// �����ת��

void vec_sign(vector<double> w, vector<double>& v);	// �����ķ��ź���

double inv_mat_infty_norm(vector<vector<double> > A);	// ��������������ƣ��Ż�����

double mat_infty_norm(vector<vector<double> > A);	// ����������

double kappa_infty(vector<vector<double> > A);	// ����������µ�������

void print_mat(vector<vector<double> > A);	// ��ӡ����

void print_vec(vector<double> b);	// ��ӡ����

double precision(vector<vector<double> > A, vector<double> b, vector<double> x);	// �����ľ��ȹ���



// =======================�������������ĺ���===================

double householder(vector<double> x, vector<double>& v);	// ����Householder�任����Ӧ���㷨3.2.1

vector<double> mat_slice_vec(vector<vector<double> > A, int i, int j, int k);	// ��������A����������A(i:j, k)

vector<vector<double> > mat_slice_mat(vector<vector<double> > A, int i1, int i2, int j1, int j2);	// ��������A�������Ӿ���A(i1:i2, j1:j2)

vector<vector<double> > mat_sub(vector<vector<double> > A, vector<vector<double> > B);	// ����������A-B����������õ��ľ���

vector<vector<double> > mat_I_sub_beta_vvT(double beta, vector<double> v);	// ����I-beta*v*v^T�����ؼ�����(����)

vector<vector<double> > mat_mul_mat(vector<vector<double> > A, vector<vector<double> > B);	// �������˷�A*B�����ؼ�����

double vec_2_norm(vector<double> x);	// ��������x��2����

vector<double> QR_decomp_householder(vector<vector<double> >& A);	// ��Householder�任ʵ�־���A��QR�ֽ⣬��Ӧ���㷨3.3.1

void QR_solver(vector<vector<double> > A, vector<double>& b);	// ����QR�ֽ⣬������Է�����Ax=b

double QR_LS(vector<vector<double> > A, vector<double> b, vector<double>& x);	// ����QR�ֽ⣬�����С�������� min ||Ax-b||_2



// =======================�������������ĺ���===================

vector<double> Jacobi_iteration(vector<double> x, vector<vector<double> > B, vector<double> g);	// Jacobi���������������1���õ�������

void Jacobi_preparation(vector<vector<double> >& A, vector<double>& b);	// ����Jacobi������ϵ������B��������g

vector<double> Gauss_Seidel_iteration(vector<double> x, vector<vector<double> > B, vector<double> g);	// Gauss-Seidel���������������1���õ�������

vector<double> SOR_iteration(vector<double> x, vector<vector<double> > A, vector<double> b, double omega);	// SOR������������1�������õ�������

vector<double> vec_sub(vector<double> a, vector<double> b);	// ������������a, b�Ĳ��a-b�����ؼ�����

double mat_flatten_2_norm(vector<vector<double> > A);	// �������AչƽΪ������2����

vector<vector<double> > pde2D_Jacobi(int n, double eps);	// Jacobi������������ָ�ʽ

double find_min(vector<vector<double> > A, int& ii, int& jj);	// �ҵ����е���С���������أ�ȡ����Сֵ���±�Ϊ(ii, jj)

vector<vector<double> > pde2D_GS(int n, double eps);	// Gauss-Seidel������������ָ�ʽ

vector<vector<double> > pde2D_SOR(int n, double eps, double omega);	// SOR������������ָ�ʽ



// =======================�������������ĺ���===================

void CGD(vector<vector<double> > A, vector<double>& x, vector<double> b, double eps, int kmax);	// �����ݶȷ�������Է�����Ax=b



// =======================�������������ĺ���===================

double powerIterationMaxEigenvalue(vector<vector<double> > A, vector<double> u);	// �ݷ������A��ģ��������ֵ�����ؼ�����

double maxAbsRoot(vector<double> coef);	// ʹ���ݷ������㲢���ض���ʽ���̵�ģ���ĸ�

vector<vector<double> > eye(int n);	// ����һ��n�׵�λ�󲢷���

void mat_put_mat(vector<vector<double> >& A, vector<vector<double> > subMat, int i, int j);	// ������subMat�������A�У�����subMat[0][0]����A[i][j]��λ��

void upperHessenberg(vector<vector<double> >& A);	// ����A����Hessenberg��

void doubleShiftQR(vector<vector<double> >& H, vector<vector<double> >& P);	// Francis˫�ز�λ�Ƶ�QR�����㷨����Ӧ�ڽ̲��㷨6.4.2

vector<vector<double> > transpose(vector<vector<double> > A);	// ���ؾ���A��ת��

bool isApproximable(vector<vector<double> > H);	// �ж�Hessenberg����H�Ƿ��Լ

bool isProposedUpTriMat(vector<vector<double> > H);	// �ж�Hessenberg����H�Ƿ�������������

void find_m_l(vector<vector<double> > H, int& m, int& l);	// �㷨6.4.3�ĸ������������(3)(ii)�Ĺ��ܣ����ҵ�����Ҫ������ķǸ�����m����С�ķǸ�����l

void implicitQR(vector<vector<double> >& A, double u);	// ����ʵ����A��ʵSchur�ֽ⣺��ʽQR�㷨����Ӧ�ڽ̲��㷨6.4.3

void implicitQREigenvalue(vector<vector<double> > A, vector<double>& Re, vector<double>& Im);	// ʹ����ʽQR�������������A��ȫ������ֵ



// =======================�������������ĺ���===================

double nonDiagNorm(vector<vector<double> > A);	// ���㷽��A�ķǶԽǡ����������䶨����̲�P211

bool isPassed(vector<vector<double> > A, double delta, int& p, int& q);	// �жϹ���Jacobi�����У���ǰ�����Ƿ񡰹��ء�

void rotateJacobi(vector<vector<double> > A, int p, int q, double& c, double& s);	// ����J(p,q,\theta)�е�\cos\theta��\sin\theta��������c, s��

vector<vector<double> > leftJacobiMul(vector<vector<double> > A, int p, int q, double c, double s);	// ��Jacobi����J(p,q,\theta)��˾���A�����ؼ�����

vector<vector<double> > rightJacobiMul(vector<vector<double> > A, int p, int q, double c, double s);	// ��Jacobi����J(p,q,\theta)�ҳ˾���A�����ؼ�����

void passingJacobiMethod(vector<vector<double> > A, vector<double>& EigenValues, vector<vector<double> >& EigenVectors, double sigma);  // ����Jacobi������Գƾ����ȫ������ֵ����������

void print_mat(vector<vector<double> > A, int digits);	// ��ӡ����A������digitsΪС�����λ�������أ�

void print_vec(vector<double> v, int digits);	// ��ӡ����v������digitsΪС�����λ�������أ�

int changeSignNum(vector<vector<double> > A, double mu);	// ���� �Գ����ԽǾ��� A��mu���ı����s_n(mu)

double bisectionMethod(vector<vector<double> > A, int m, double eps);	// ʹ�ö��ַ��������A�ĵ�mС������ֵ

vector<double> inversePowerMethod(vector<vector<double> > A, double lambda, double eps);	// ʹ�÷��ݷ�������A������ֵlambda��Ӧ����������

void helper(vector<double> eigen);	// ����С����



// =======================�ڰ����������ĺ���===================

vector<vector<double> > xyT(vector<double> x, vector<double> y);	// ������� xy^T������x, yΪ����������С�ֱ�Ϊm, n

vector<double> numMul(double lambda, vector<double> v);	// �����������ˣ����ؼ�����

vector<vector<double> > numMul(double lambda, vector<vector<double> > A);	// ����������ˣ����ؼ�����

vector<vector<double> > biDiagonalization(vector<vector<double> > A, vector<vector<double> >& U, vector<vector<double> >& V);	// �Ծ��� A �� Householder �任�����Խǻ�

void givens(double& c, double& s, double a, double b);	// ��Givens�任�е�\cos\theta��\sin\theta���ֱ𱣴���c, s��

vector<vector<double> > WilkinsonSVD(vector<vector<double> > B, vector<vector<double> >& P, vector<vector<double> >& Q);	// ��Wilkinsonλ�Ƶ�SVD����

void find_p_q(vector<vector<double> > B, int& p, int& q);	// �㷨7.6.3�еĸ�������������Ҫ���p��q

void SVD(vector<vector<double> > A, vector<vector<double> >& U, vector<vector<double> >& Sigma, vector<vector<double> >& V, double eps);	// ����A��SVD�ֽ⣬����A�Ĵ�СΪm*n��m\geq n

double maxAbs(vector<vector<double> > A);	// ���ؾ���A��ģ����Ԫ��