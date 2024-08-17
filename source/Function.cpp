#include "Function.h"
#include <cmath>
#include <iomanip>
#include <algorithm>

// ============================第一章中定义的函数=============================

void forward_subs(vector<vector<double> > L, vector<double>& b)
{
    int n = L.size();
    for (int j = 0; j < n - 1; j++) {
        b[j] /= L[j][j];
        for (int i = j + 1; i < n; i++) {
            b[i] -= b[j] * L[i][j];
        }
    }
    b[n - 1] /= L[n - 1][n - 1];
}

void forward_subs1(vector<vector<double> > L, vector<double>& b)
{
    int n = L.size();
    for (int j = 0; j < n - 1; j++) {
        for (int i = j + 1; i < n; i++) {
            b[i] -= b[j] * L[i][j];
        }
    }
}

void back_subs(vector<vector<double> > U, vector<double>& y)
{
    int n = U[0].size();
    for (int j = n - 1; j > 0; j--) {
        y[j] /= U[j][j];
        for (int i = 0; i < j; i++)
            y[i] -= y[j] * U[i][j];
    }
    y[0] /= U[0][0];
}

void back_subs1(vector<vector<double> > U, vector<double>& y)
{
    int n = U.size();
    for (int j = n - 1; j > 0; j--) {
        for (int i = 0; i < j; i++)
            y[i] -= y[j] * U[i][j];
    }
}

void gauss_elim(vector<vector<double> >& A)
{
    int n = A.size();
    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {
            A[i][k] /= A[k][k];
        }
        for (int i = k + 1; i < n; i++) {
            for (int j = k + 1; j < n; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
}

void gauss_elim_full_pivoting(vector<vector<double> >& A, vector<int>& u, vector<int>& v)
{
    int n = A.size();
    int p, q;
    double temp_max;
    for (int k = 0; k < n - 1; k++) {

        // find the argmax p, q
        p = k; q = k; temp_max = fabs(A[p][q]);
        for (int i = k; i < n; i++) {
            for (int j = k; j < n; j++) {
                if (fabs(A[i][j]) > temp_max) {
                    p = i;
                    q = j;
                    temp_max = fabs(A[p][q]);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            std::swap(A[k][i], A[p][i]);
        }
        for (int i = 0; i < n; i++) {
            std::swap(A[i][k], A[i][q]);
        }

        u.push_back(p);
        v.push_back(q);

        if (A[k][k] != 0) {
            for (int i = k + 1; i < n; i++) {
                A[i][k] /= A[k][k];
            }
            for (int i = k + 1; i < n; i++) {
                for (int j = k + 1; j < n; j++) {
                    A[i][j] -= A[i][k] * A[k][j];
                }
            }
        }
        else {
            cout << "The matrix A is singular!" << endl;
            return;
        }
    }
    u.push_back(n - 1);
    v.push_back(n - 1);
}

void gauss_elim_col_pivoting(vector<vector<double> >& A, vector<int>& u)
{
    int n = A.size();
    int p;
    double temp_max;
    u.clear();
    for (int k = 0; k < n - 1; k++) {
        p = k;
        temp_max = fabs(A[k][k]);
        for (int i = k; i < n; i++) {
            if (fabs(A[i][k]) > temp_max) {
                temp_max = fabs(A[i][k]);
                p = i;
            }
        }
        for (int i = 0; i < n; i++) {
            std::swap(A[k][i], A[p][i]);
        }
        u.push_back(p);
        if (A[k][k]) {
            for (int i = k + 1; i < n; i++) {
                A[i][k] /= A[k][k];
            }
            for (int i = k + 1; i < n; i++) {
                for (int j = k + 1; j < n; j++) {
                    A[i][j] -= A[i][k] * A[k][j];
                }
            }
        }
        else {
            cout << "The matrix A is singular!" << endl;
            return;
        }
    }
    u.push_back(n - 1);
}

void vector_pb(vector<int>& u, vector<double>& b)
{
    int n = u.size();
    for (int i = 0; i < n; i++) {
        std::swap(b[i], b[u[i]]);
    }
}

void vector_qb(vector<int>& v, vector<double>& b)
{
    int n = v.size();
    for (int i = 0; i < n; i++) {
        std::swap(b[i], b[v[i]]);
    }
}

void cholesky_decomp(vector<vector<double> >& A)
{
    int n = A.size();
    for (int k = 0; k < n; k++) {
        A[k][k] = sqrt(A[k][k]);
        for (int i = k + 1; i < n; i++) {
            A[i][k] /= A[k][k];
        }
        for (int j = k + 1; j < n; j++) {
            for (int i = j; i < n; i++) {
                A[i][j] -= A[i][k] * A[j][k];
            }
        }
    }
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            A[i][j] = A[j][i];
        }
    }
}

void modified_cholesky_decomp(vector<vector<double> >& A)
{
    int n = A.size();
    vector<double> v;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < j; i++) {
            v.push_back(A[j][i] * A[i][i]);
        }
        for (int k = 0; k < j; k++) {
            A[j][j] -= A[j][k] * v[k];
        }
        for (int k = j + 1; k < n; k++) {
            for (int i = 0; i < j; i++) {
                A[k][j] -= A[k][i] * v[i];
            }
            A[k][j] /= A[j][j];
        }
        v.clear();
    }
}

void matrix_DLT(vector<vector<double> >& A)
{
    int n = A.size();
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            A[i][j] = A[j][i] * A[i][i];
        }
    }
}



// =============================第二章中定义的函数=============================

void mat_mul_vec(vector<vector<double> > A, vector<double> x, vector<double>& b) {
    // 计算矩阵A和向量x的乘积，计算结果保存在向量b中
    // 向量b要求是大小和x相同的vector<double>
    int n = x.size();
    double sum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        b[i] = sum;
        sum = 0;
    }
}

void sol_eqs_col_LU(vector<vector<double> > A, vector<double> b, vector<double>& sol) {//checked
    // 使用列主元Gauss消去法求解线性方程组Ax=b，将计算结果保存在向量sol中
    // 传入的参数sol指定为大小和b相同的vector<double>
    int n = b.size();
    vector<int> u;
    gauss_elim_col_pivoting(A, u);
    vector_pb(u, b);
    forward_subs1(A, b);
    back_subs(A, b);
    for (int i = 0; i < n; i++) {
        sol[i] = b[i];
    }
}

double vec_1_norm(vector<double> w) {//checked
    // 计算向量w的1范数
    double norm = 0;
    for (vector<double>::iterator it = w.begin(); it != w.end(); it++) {
        norm += fabs(*it);
    }
    return norm;
}

double vec_infty_norm(vector<double> z, int &j) {//checked
    // 计算向量z的无穷范数和使得|z_j|=|z|_\infty的下标j
    double norm = 0; j = 0;
    for (int i = 0; i < z.size(); i++) {
        if (fabs(z[i]) > norm){
            norm = fabs(z[i]);
            j = i;
        }
    }
    return norm;
}

double vec_infty_norm(vector<double> z) {
    // 计算向量z的无穷范数（重载）
    double norm = 0;
    for (int i = 0; i < z.size(); i++) {
        if (fabs(z[i]) > norm) {
            norm = fabs(z[i]);
        }
    }
    return norm;
}

double vec_inner_prod(vector<double> x, vector<double> y) {//checked
    // 计算向量x, y的内积
    int n = x.size();
    double inner_prod = 0;
    for (int i = 0; i < n; i++) {
        inner_prod += x[i] * y[i];
    }
    return inner_prod;
}

void transpose(vector<vector<double> > A, vector<vector<double> >& At) {//checked
    // 计算矩阵A的转置，并存入At中
    int n = A.size();
    int m = A[0].size();
    // A的大小为n*m, A^T的大小为m*n
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            At[j][i] = A[i][j];
        }
    }
}

void vec_sign(vector<double> w, vector<double> &v) {//checked
    // 计算sign(w)，结果存入v中
    // v要求是大小和w相同的vector<double>
    int n = w.size();
    for (int i = 0; i < n; i++) {
        v[i] = ((w[i]>0) - (w[i]<0));
    }
}

double inv_mat_infty_norm(vector<vector<double> > A) {
    // 算法2.5.1的修改：估计矩阵A的逆矩阵的无穷范数
    int k = 1; 
    int n = A.size();
    vector<double> x;

    for (int i = 0; i < n; i++) {   // 初始化向量x
        x.push_back(1.0 / n);
    }

    vector<double> w(n);
    vector<double> z(n);
    vector<double> v(n);
    vector<vector<double> > At(n, vector<double>(n));
    transpose(A, At);
    int j = 0;
    double norm;
    
    while (k) {
        sol_eqs_col_LU(At, x, w);
        vec_sign(w, v);
        sol_eqs_col_LU(A, w, z);

        norm = vec_infty_norm(z, j);

        if (norm <= vec_inner_prod(z, x)) {
            return vec_1_norm(w);
        }
        else {
            for (int i = 0; i < n; i++) {
                x[i] = 0;
            }
            x[j] = 1;
            k += 1;
        }
    }
    return -1;
}

double mat_infty_norm(vector<vector<double> > A) {//checked
    // 计算矩阵A的无穷范数
    int n = A.size();
    double norm = 0;
    double row = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            row += fabs(A[i][j]);
        }
        if (norm < row) {
            norm = row;
        }
        row = 0;
    }
    return norm;
}

double kappa_infty(vector<vector<double> > A) {
    // 计算矩阵A的无穷范数意义下的条件数    
    return inv_mat_infty_norm(A) * mat_infty_norm(A);
}

void print_mat(vector<vector<double> > A) {
    // 打印矩阵A
    int m = A.size();   // A为m行
    int n = A[0].size();// A为n列
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << '\t';
        }
        cout << endl;
    }
}

void print_vec(vector<double> b) {
    // 打印向量b，以行向量的形式打印
    for (vector<double>::iterator it = b.begin(); it != b.end(); it++) {
        cout << *it << '\t';
    }
    cout << endl;
}

double precision(vector<vector<double> > A, vector<double> b, vector<double> x) {//checked
    // 在无穷范数意义下，估计计算解的精度
    // A是系数矩阵，b是右端项，x是计算解
    int n = b.size();
    vector<double> r(n);
    mat_mul_vec(A, x, r);
    
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i]; 
    }

    int idx;
    double r_infty = vec_infty_norm(r, idx);
    double b_infty = vec_infty_norm(b, idx);

    return kappa_infty(A) * r_infty / b_infty;
}



// =============================第三章中定义的函数=============================

double householder(vector<double> x, vector<double>& v) {
    // 计算Householder变换，对应于算法3.2.1
    // 输入x为原来的向量，v用于存储Householder变换中的向量v(需要事先指定大小)，返回值为beta
    int n = x.size();
    int idx = 0;
    double eta = vec_infty_norm(x, idx);
    double beta;
    for (int i = 0; i < n; i++) {
        x[i] /= eta;
    }
    double sigma = 0;
    for (int i = 1; i < n; i++) {
        sigma += x[i] * x[i];
        v[i] = x[i];
    }
    if (!sigma) {
        beta = 0;
    }
    else {
        double alpha = sqrt(x[0] * x[0] + sigma);
        if (x[0] <= 0) {
            v[0] = x[0] - alpha;
        }
        else {
            v[0] = -sigma / (x[0] + alpha);
        }
        beta = 2 * v[0] * v[0] / (sigma + v[0] * v[0]);
        double d = v[0];
        for (int i = 0; i < n; i++) {
            v[i] /= d;
        }
    }
    return beta;
}

vector<double> mat_slice_vec(vector<vector<double> > A, int i, int j, int k) {
    // 给定矩阵A，返回向量A(i:j, k)
    vector<double> result;
    for (int p = i; p < j + 1; p++) {
        result.push_back(A[p][k]);
    }
    return result;
}

vector<vector<double> > mat_slice_mat(vector<vector<double> > A, int i1, int i2, int j1, int j2) {
    // 给定矩阵A，返回子矩阵A(i1:i2, j1:j2)
    int row = i2 - i1 + 1;  // 子矩阵的行数
    int col = j2 - j1 + 1;  // 子矩阵的列数
    vector<vector<double> > result(row, vector<double>(col));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            result[i][j] = A[i1 + i][j1 + j];
        }
    }
    return result;
}

vector<vector<double> > mat_sub(vector<vector<double> > A, vector<vector<double> > B) {
    // 计算矩阵减法A-B，返回相减得到的矩阵
    // 要求A，B相同大小
    int m = A.size();       // 行数
    int n = A[0].size();    // 列数
    vector<vector<double> > result(m, vector<double>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

vector<vector<double> > mat_I_sub_beta_vvT(double beta, vector<double> v) {
    // 计算I-beta*v*v^T，返回计算结果(矩阵)
    int n = v.size();
    vector<vector<double> > A(n, vector<double>(n));    // 用于存储I
    vector<vector<double> > B(n, vector<double>(n));    // 用于存储beta*v*v^T
    for (int i = 0; i < n; i++) {
        A[i][i] = 1;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B[i][j] = beta * v[i] * v[j];
        }
    }
    vector<vector<double> > result = mat_sub(A, B);
    return result;
}

vector<vector<double> > mat_mul_mat(vector<vector<double> > A, vector<vector<double> > B) {
    // 计算矩阵乘法A*B，返回计算结果
    // 这里要求乘法是相容的，即A大小为p*q，B大小为q*r，计算结果的大小为p*r
    int p = A.size();
    int q = A[0].size();
    int r = B[0].size();
    vector<vector<double> > result(p, vector<double>(r));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < r; j++) {
            for (int k = 0; k < q; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

double vec_2_norm(vector<double> x) {
    // 计算向量x的2范数
    int n = x.size();
    double norm = 0;
    for (int i = 0; i < n; i++) {
        norm += x[i] * x[i];
    }
    norm = sqrt(norm);
    return norm;
}

vector<double> QR_decomp_householder(vector<vector<double> >& A) {
    // 用Householder变换实现矩阵A的QR分解，对应于算法3.3.1
    // 分解得到的Q用n个Householder变换的v和beta表示，存储在A的上三角元以外的元素中；R存储在A的上三角元中
    int m = A.size();       // 矩阵A有m行
    int n = A[0].size();    // 矩阵A有n列
    vector<double> d;
    for (int j = 0; j < n; j++) {
        if (j < m - 1) {
            vector<double> x = mat_slice_vec(A, j, m - 1, j);
            vector<double> v(m-j);  // v的大小是m-j
            double beta = householder(x, v);
            vector<vector<double> > temp1 = mat_slice_mat(A, j, m - 1, j, n - 1); // temp1的大小是(m-j)*(n-j)
            vector<vector<double> > temp2 = mat_I_sub_beta_vvT(beta, v);    // temp2 = I-beta*v*v^T
            temp1 = mat_mul_mat(temp2, temp1);  // temp1 = (I-beta*v*v^T)*A(j:m-1,j:n-1)
            for (int p = 0; p < m - j; p++) {
                for (int q = 0; q < n - j; q++) {
                    A[p + j][q + j] = temp1[p][q];
                }
            }
            d.push_back(beta);
            for (int p = j + 1; p < m; p++) {
                A[p][j] = v[p - j];
            }
        }
    }
    return d;
}

void QR_solver(vector<vector<double> > A, vector<double>& b) {
    // 利用QR分解，求解线性方程组Ax=b
    // 要求A是可逆方阵，求出的解保存在向量b中
    int n = A.size();
    vector<double> d = QR_decomp_householder(A);
    vector<double> v(n);
    vector<vector<double> > H(n, vector<double>(n));
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < i; j++) {
            v[j] = 0;
        }
        v[i] = 1;
        for (int j = i + 1; j < n; j++) {
            v[j] = A[j][i];
        }
        H = mat_I_sub_beta_vvT(d[i], v);
        mat_mul_vec(H, b, b);
    }
    back_subs(A, b);
}

double QR_LS(vector<vector<double> > A, vector<double> b, vector<double>& x) {
    // 利用QR分解，求解最小二乘问题 min ||Ax-b||_2
    // 返回上述最小二乘的值，最小二乘解x保存在向量x中，其中传入的x的大小等于A的列数；要求A的行数严格大于列数
    int m = A.size();       // 矩阵A的行数
    int n = A[0].size();    // 矩阵A的列数
    vector<double> d = QR_decomp_householder(A);
    vector<double> v(m);
    vector<vector<double> > H(m, vector<double>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            v[j] = 0;
        }
        v[i] = 1;
        for (int j = i + 1; j < m; j++) {
            v[j] = A[j][i];
        }
        H = mat_I_sub_beta_vvT(d[i], v);
        mat_mul_vec(H, b, b);
    }
    vector<double> c1(n), c2(m - n);
    for (int i = 0; i < n; i++) {
        c1[i] = b[i];
    }
    for (int i = n; i < m; i++) {
        c2[i - n] = b[i];
    }
    back_subs(A, c1);
    for (int i = 0; i < c1.size(); i++) {
        x[i] = c1[i];
    }
    return vec_2_norm(c2);
}



// =============================第四章中定义的函数=============================

vector<double> Jacobi_iteration(vector<double> x, vector<vector<double> > B, vector<double> g) {
    // Jacobi迭代法，计算迭代1步得到的向量，其中 x_{k+1} = Bx_k + g
    // x-迭代前的向量x；B-迭代矩阵；g-常数向量；返回-迭代1步后的向量y
    int n = x.size();
    vector<double> y(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            y[i] += B[i][j] * x[j];
        }
        y[i] += g[i];
    }
    return y;
}

void Jacobi_preparation(vector<vector<double> >& A, vector<double>& b) {
    // 计算Jacobi迭代的系数矩阵B、常向量g
    // 计算结果分别保存在原始的系数矩阵A、左端向量b中
    int n = b.size();
    vector<vector<double> > B = A;
    for (int i = 0; i < n; i++) {   // 计算常数向量g
        b[i] /= A[i][i];
    }
    for (int i = 0; i < n; i++) {   // 计算迭代矩阵B
        for (int j = 0; j < n; j++) {
            B[i][j] = -A[i][j] / A[i][i];
        }
    }
    for (int i = 0; i < n; i++) {   // 修正矩阵B的对角元为0
        B[i][i] = 0;
    }
    A = B;
}

vector<double> Gauss_Seidel_iteration(vector<double> x, vector<vector<double> > B, vector<double> g) {
    // Gauss-Seidel迭代法，计算迭代1步得到的向量
    // x-迭代前的向量x；B-迭代矩阵；g-常数向量；返回-迭代1步后的向量y
    int n = x.size();
    for (int i = 0; i < n; i++) {
        x[i] = 0;
        for (int j = 0; j < n; j++) {
            x[i] += B[i][j] * x[j];
        }
        x[i] += g[i];
    }
    return x;
}

vector<double> SOR_iteration(vector<double> x, vector<vector<double> > A, vector<double> b, double omega) {
    // SOR迭代法，计算1步迭代得到的向量
    // x-迭代前的向量x；A-系数矩阵；b-右端向量；omega-松弛因子
    int n = x.size();
    vector<double> c(n);
    vector<vector<double> > B(n, vector<double>(n));    // 先存放(1-omega)D + omega U，之后存放D - omega L
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {   // 计算矩阵B
            B[i][j] = -omega * A[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        B[i][i] = (1 - omega) * A[i][i];
    }   // 现在B存储的是(1-omega)D + omega U
    mat_mul_vec(B, x, c);   // 计算矩阵向量乘法：c = Bx
    for (int i = 0; i < n; i++) {
        b[i] = omega * b[i] + c[i];
    }   // 现在b存储的是[(1-omega)D+omega U]x + omega b
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            B[i][j] = 0;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = i - 1; j >= 0; j--) {
            B[i][j] = omega * A[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        B[i][i] = A[i][i];
    }   // 现在B = (D - omega L)
    forward_subs(B, b);
    return b;
}

vector<double> vec_sub(vector<double> a, vector<double> b) {
    // 计算两个向量a, b的差，即a-b，返回计算结果
    // 这里要求这两个向量大小相同
    vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

double mat_flatten_2_norm(vector<vector<double> > A) {
    // 计算矩阵A展平为向量的2范数
    int m = A.size();
    int n = A[0].size();
    double norm = 0.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            norm += A[i][j] * A[i][j];
        }
    }
    return sqrt(norm);
}

vector<vector<double> > pde2D_Jacobi(int n, double eps) {
    // Jacobi迭代法，五点差分格式，求解pde: -\nabla^2 u + g(x,y)u = f(x,y), (x,y)\in [0,1]^2; 边值 u=1
    // 这里 g(x,y) = exp(xy), f(x,y) = x + y
    // n-网格线数目；eps-停机条件
    // 返回求出的数值解
    double h = 1.0 / n;
    vector<vector<double> > uCur(n - 1, vector<double>(n - 1)), uPre(n - 1, vector<double>(n - 1, 1));
    // uCur是迭代1步后的解，uPre是迭代1步前的解
    int iterNum = 0;    // 迭代次数
    
    double a, b;
    while (mat_flatten_2_norm(mat_sub(uCur, uPre)) > eps) {
        uPre = uCur;
        for (int i = 1; i < n - 2; i++) {
            a = 4 + h * h * exp(h * h * (i + 1));
            b = (i + 2) * h * h * h;
            uCur[i][0] = (1 + uPre[i][1] + uPre[i + 1][0] + uPre[i - 1][0] + b) / a;
        }
        for (int i = 1; i < n - 2; i++) {
            a = 4 + h * h * exp((n - 1) * (i + 1) * h * h);
            b = (i + n) * h * h * h;
            uCur[i][n - 2] = (1 + uPre[i][n - 3] + uPre[i - 1][n - 2] + uPre[i + 1][n - 2] + b) / a;
        }
        for (int j = 1; j < n - 2; j++) {
            a = 4 + h * h * exp((j + 1) * h * h);
            b = (j + 2) * h * h * h;
            uCur[0][j] = (1 + uPre[1][j] + uPre[0][j - 1] + uPre[0][j + 1] + b) / a;
        }
        for (int j = 1; j < n - 2; j++) {
            a = 4 + h * h * exp((n - 1) * (j + 1) * h * h);
            b = (n + j) * h * h * h;
            uCur[n - 2][j] = (1 + uPre[n - 3][j] + uPre[n - 2][j - 1] + uPre[n - 2][j + 1] + b) / a;
        }
        
        a = 4 + h * h * exp(h * h);
        b = 2 * h * h * h;
        uCur[0][0] = (1 + 1 + uPre[0][1] + uPre[1][0] + b) / a;
        
        a = 4 + h * h * exp((n - 1) * h * h);
        b = n * h * h * h;
        uCur[n - 2][0] = (1 + 1 + uPre[n - 2][1] + uPre[n - 3][0] + b) / a;
        
        a = 4 + h * h * exp((n - 1) * (n - 1) * h * h);
        b = 2 * (n - 1) * h * h * h;
        uCur[n - 2][n - 2] = (1 + 1 + uPre[n - 2][n - 3] + uPre[n - 3][n - 2] + b) / a;
        
        a = 4 + h * h * exp((n - 1) * h * h);
        b = n * h * h * h;
        uCur[0][n - 2] = (1 + 1 + uPre[0][n - 3] + uPre[1][n - 2] + b) / a;

        for (int i = 1; i < n - 2; i++) {
            for (int j = 1; j < n - 2; j++) {
                a = 4 + h * h * exp((i + 1) * (j + 1) * h * h);
                b = h * h * (i + j + 2) * h;
                uCur[i][j] = uPre[i - 1][j] / a + uPre[i][j - 1] / a 
                    + uPre[i + 1][j] / a + uPre[i][j + 1] / a + b / a;
            }
        }

        iterNum++;
    }

    cout << "iteration times: " << iterNum << endl;
    return uCur;
}

double find_min(vector<vector<double> > A, int& ii, int& jj) {
    // 找到解中的最小分量并返回，取到最小值的下标为(ii, jj)
    // 默认A是方阵
    int n = A.size();
    double min = INFINITY;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A[i][j] < min) {
                min = A[i][j];
                ii = i;
                jj = j;
            }
        }
    }
    return min;
}

vector<vector<double> > pde2D_GS(int n, double eps) {
    // Gauss-Seidel迭代法，五点差分格式，求解pde: -\nabla^2 u + g(x,y)u = f(x,y), (x,y)\in [0,1]^2; 边值 u=1
    // 这里 g(x,y) = exp(xy), f(x,y) = x + y
    // n-网格线数目；eps-停机条件
    // 返回求出的数值解
    double h = 1.0 / n;
    vector<vector<double> > uCur(n - 1, vector<double>(n - 1)), uPre(n - 1, vector<double>(n - 1, 1));
    // uCur是迭代1步后的解，uPre是迭代1步前的解
    int iterNum = 0;    // 迭代次数

    double a, b;
    while (mat_flatten_2_norm(mat_sub(uCur, uPre)) > eps) {

        uPre = uCur;
        
        a = 4 + h * h * exp(h * h);
        b = 2 * h * h * h;
        uCur[0][0] = (1 + 1 + uCur[0][1] + uCur[1][0] + b) / a;

        for (int i = 1; i < n - 2; i++) {
            a = 4 + h * h * exp(h * h * (i + 1));
            b = (i + 2) * h * h * h;
            uCur[i][0] = (1 + uCur[i][1] + uCur[i + 1][0] + uCur[i - 1][0] + b) / a;
        }

        a = 4 + h * h * exp((n - 1) * h * h);
        b = n * h * h * h;
        uCur[n - 2][0] = (1 + 1 + uCur[n - 2][1] + uCur[n - 3][0] + b) / a;
        

        for (int j = 1; j < n - 2; j++) {
            a = 4 + h * h * exp((j + 1) * h * h);
            b = (j + 2) * h * h * h;
            uCur[0][j] = (1 + uCur[1][j] + uCur[0][j - 1] + uCur[0][j + 1] + b) / a;
            
            for (int i = 1; i < n - 2; i++) {
                a = 4 + h * h * exp((i + 1) * (j + 1) * h * h);
                b = h * h * (i + j + 2) * h;
                uCur[i][j] = (uCur[i - 1][j] + uCur[i][j - 1] 
                    + uCur[i + 1][j] + uCur[i][j + 1] + b) / a;
            }
            
            a = 4 + h * h * exp((n - 1) * (j + 1) * h * h);
            b = (n + j) * h * h * h;
            uCur[n - 2][j] = (1 + uCur[n - 3][j] + uCur[n - 2][j - 1] + uCur[n - 2][j + 1] + b) / a;
        }


        a = 4 + h * h * exp((n - 1) * h * h);
        b = n * h * h * h;
        uCur[0][n - 2] = (1 + 1 + uCur[0][n - 3] + uCur[1][n - 2] + b) / a;

        for (int i = 1; i < n - 2; i++) {
            a = 4 + h * h * exp((n - 1) * (i + 1) * h * h);
            b = (i + n) * h * h * h;
            uCur[i][n - 2] = (1 + uCur[i][n - 3] + uCur[i - 1][n - 2] + uCur[i + 1][n - 2] + b) / a;
        }
        
        a = 4 + h * h * exp((n - 1) * (n - 1) * h * h);
        b = 2 * (n - 1) * h * h * h;
        uCur[n - 2][n - 2] = (1 + 1 + uCur[n - 2][n - 3] + uCur[n - 3][n - 2] + b) / a;

        iterNum++;
    }

    cout << "iteration times: " << iterNum << endl;
    return uCur;
}

vector<vector<double> > pde2D_SOR(int n, double eps, double omega) {
    // SOR迭代法，五点差分格式，求解pde: -\nabla^2 u + g(x,y)u = f(x,y), (x,y)\in [0,1]^2; 边值 u=1
    // 这里 g(x,y) = exp(xy), f(x,y) = x + y
    // n-网格线数目；eps-停机条件；omega-松弛因子
    // 返回求出的数值解
    double h = 1.0 / n;
    vector<vector<double> > uCur(n - 1, vector<double>(n - 1)), uPre(n - 1, vector<double>(n - 1, 1));
    // uCur是迭代1步后的解，uPre是迭代1步前的解
    int iterNum = 0;    // 迭代次数

    double a, b;
    while (mat_flatten_2_norm(mat_sub(uCur, uPre)) > eps) {

        uPre = uCur;

        a = 4 + h * h * exp(h * h);
        b = 2 * h * h * h;
        uCur[0][0] = (1 - omega) * uCur[0][0] + omega * (1 + 1 + uCur[0][1] + uCur[1][0] + b) / a;

        for (int i = 1; i < n - 2; i++) {
            a = 4 + h * h * exp(h * h * (i + 1));
            b = (i + 2) * h * h * h;
            uCur[i][0] = (1 - omega) * uCur[i][0]
                + omega * (1 + uCur[i][1] + uCur[i + 1][0] + uCur[i - 1][0] + b) / a;
        }

        a = 4 + h * h * exp((n - 1) * h * h);
        b = n * h * h * h;
        uCur[n - 2][0] = (1 - omega) * uCur[n - 2][0]
            + omega * (1 + 1 + uCur[n - 2][1] + uCur[n - 3][0] + b) / a;


        for (int j = 1; j < n - 2; j++) {
            a = 4 + h * h * exp((j + 1) * h * h);
            b = (j + 2) * h * h * h;
            uCur[0][j] = (1 - omega) * uCur[0][j]
                + omega * (1 + uCur[1][j] + uCur[0][j - 1] + uCur[0][j + 1] + b) / a;

            for (int i = 1; i < n - 2; i++) {
                a = 4 + h * h * exp((i + 1) * (j + 1) * h * h);
                b = h * h * (i + j + 2) * h;
                uCur[i][j] = (1 - omega) * uCur[i][j]
                    + omega * (uCur[i - 1][j] + uCur[i][j - 1] + uCur[i + 1][j] + uCur[i][j + 1] + b) / a;
            }

            a = 4 + h * h * exp((n - 1) * (j + 1) * h * h);
            b = (n + j) * h * h * h;
            uCur[n - 2][j] = (1 - omega) * uCur[n - 2][j]
                + omega * (1 + uCur[n - 3][j] + uCur[n - 2][j - 1] + uCur[n - 2][j + 1] + b) / a;
        }


        a = 4 + h * h * exp((n - 1) * h * h);
        b = n * h * h * h;
        uCur[0][n - 2] = (1 - omega) * uCur[0][n - 2]
            + omega * (1 + 1 + uCur[0][n - 3] + uCur[1][n - 2] + b) / a;

        for (int i = 1; i < n - 2; i++) {
            a = 4 + h * h * exp((n - 1) * (i + 1) * h * h);
            b = (i + n) * h * h * h;
            uCur[i][n - 2] = (1 - omega) * uCur[i][n - 2]
                + omega * (1 + uCur[i][n - 3] + uCur[i - 1][n - 2] + uCur[i + 1][n - 2] + b) / a;
        }

        a = 4 + h * h * exp((n - 1) * (n - 1) * h * h);
        b = 2 * (n - 1) * h * h * h;
        uCur[n - 2][n - 2] = (1 - omega) * uCur[n - 2][n - 2]
            + omega * (1 + 1 + uCur[n - 2][n - 3] + uCur[n - 3][n - 2] + b) / a;

        iterNum++;
    }

    cout << "iteration times: " << iterNum << endl;
    return uCur;
}



// =============================第五章中定义的函数=============================

void CGD(vector<vector<double> > A, vector<double>& x, vector<double> b, double eps, int kmax) {
    // 共轭梯度法求解线性方程组Ax=b
    // A-系数矩阵；x-初值；b-右端项；eps-停机条件；kmax-最大迭代次数；最终求出的解存储在x中
    int k = 0;
    vector<double> r(x.size());
    vector<double> p(x.size());
    mat_mul_vec(A, x, r);
    r = vec_sub(b, r);
    double rho = vec_inner_prod(r, r);
    double stop = pow(eps * vec_2_norm(b), 2);
    double beta, rhoTilde, alpha;
    vector<double> w(x.size());
    while ((rho > stop) && (k < kmax)) {
        k++;
        if (k == 1) {
            p = r;
        }
        else {
            beta = rho / rhoTilde;
            for (int i = 0; i < x.size(); i++) {
                p[i] = r[i] + beta * p[i];
            }
        }
        mat_mul_vec(A, p, w);
        alpha = rho / vec_inner_prod(p, w);
        for (int i = 0; i < x.size(); i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * w[i];
        }
        rhoTilde = rho;
        rho = vec_inner_prod(r, r);
    }
    cout << "iteration: " << k << endl;
}



// =============================第六章中定义的函数=============================

double powerIterationMaxEigenvalue(vector<vector<double> > A, vector<double> u) {
    // 幂法求矩阵A的模最大的特征值，返回计算结果
    // A-待求特征值的矩阵；u-迭代初始向量，要求||u||_\infty = 1
    vector<double> y(u.size());
    double mu = 1;
    double muPre = 0;
    int idx;
    int count = 0;
    int maxIteration = 1000;    // 最大迭代次数：1000
    while (abs(mu - muPre) > 1e-7) {  // 迭代终止容差：1e-7
        muPre = mu;
        mat_mul_vec(A, u, y);
        vec_infty_norm(y, idx);
        mu = y[idx];
        for (int i = 0; i < u.size(); i++) {
            u[i] = y[i] / mu;
        }
        count++;
        if (count > maxIteration) {
            cout << "Max iteration reached!" << endl;
            break;
        }
    }
    cout << "Power Iterations: " << count << endl;
    return mu;
}

double maxAbsRoot(vector<double> coef) {
    // 使用幂法，计算并返回多项式方程的模最大的根
    // 多项式方程为：x^n + coef[n-1]x^{n-1} + ... + coef[1]x + coef[0] = 0
    // coef: 多项式的系数
    int n = coef.size();    // 多项式的次数
    vector<vector<double> > A(n, vector<double>(n));    // 多项式的友阵
    vector<double> u(n); // 迭代的初始向量
    for (int i = 0; i < n; i++) {
        A[0][i] = -coef[n - 1 - i];
    }
    for (int i = 1; i < n; i++) {
        A[i][i - 1] = 1;
    }
    u[0] = -0.9;
    u[n - 1] = 0.9;
    double root = powerIterationMaxEigenvalue(A, u);
    return root;
}

vector<vector<double> > eye(int n) {
    // 生成一个n阶单位阵并返回
    vector<vector<double> > I(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        I[i][i] = 1;
    }
    return I;
}

void mat_put_mat(vector<vector<double> >& A, vector<vector<double> > subMat, int i, int j) {
    // 将矩阵subMat放入矩阵A中，其中subMat[0][0]放在A[i][j]的位置
    if ((subMat.size() + i > A.size()) || (subMat[0].size() + j > A[0].size())) {
        cout << "ERROR: Matrix index out of range!" << endl;
        return;
    }
    for (int x = 0; x < subMat.size(); x++) {
        for (int y = 0; y < subMat[0].size(); y++) {
            A[x + i][y + j] = subMat[x][y];
        }
    }
}

void upperHessenberg(vector<vector<double> >& A) {
    // 矩阵A的上Hessenberg化，对应于教材算法6.4.1
    // 得到A的上Hessenberg矩阵仍存储在A中
    int n = A.size();
    vector<double> v;
    vector<double> x;
    vector<vector<double> > subMat;
    double beta;
    int idx;
    for (int k = 0; k < n - 2; k++) {
        x = mat_slice_vec(A, k + 1, n - 1, k);
        if (vec_infty_norm(x, idx) == 0) {
            continue;
        }
        v.resize(n - k - 1);
        beta = householder(x, v);
        subMat = mat_mul_mat(mat_I_sub_beta_vvT(beta, v), mat_slice_mat(A, k+1, n-1, k, n-1));
        mat_put_mat(A, subMat, k + 1, k);
        subMat.clear();

        subMat = mat_mul_mat(mat_slice_mat(A, 0, n-1, k+1, n-1), mat_I_sub_beta_vvT(beta, v));
        mat_put_mat(A, subMat, 0, k + 1);

        x.clear();
        v.clear();
        subMat.clear();
    }

    for (int i = 2; i < n; i++) {
        for (int j = 0; j < i - 1; j++) {
            A[i][j] = 0;
        }
    }

}

void doubleShiftQR(vector<vector<double> >& H, vector<vector<double> >& P) {
    // Francis双重步位移的QR迭代算法，对应于教材算法6.4.2
    // 这里H要求是一个上Hessenberg矩阵，计算迭代1步后的结果，并保存在H中
    // P用于存储n-1个Householder变换的乘积
    int n = H.size();
    if (n > 2) {
        int m = n - 1;
        double s = H[m - 1][m - 1] + H[n - 1][n - 1];
        double t = H[m - 1][m - 1] * H[n - 1][n - 1] - H[m - 1][n - 1] * H[n - 1][m - 1];
        double x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
        double y = H[1][0] * (H[0][0] + H[1][1] - s);
        double z = H[1][0] * H[2][1];
        vector<double> v(3);
        vector<double> u(3);
        double beta;
        int q;
        int r;
        vector<vector<double> > subMat;
        P = eye(n); // 初始化 P = I_n
        for (int k = 0; k < n - 2; k++) {
            u[0] = x; u[1] = y; u[2] = z;
            beta = householder(u, v);
            
            q = max(1, k);
            subMat = mat_mul_mat(mat_I_sub_beta_vvT(beta, v), mat_slice_mat(H, k, k + 2, q - 1, n - 1));
            mat_put_mat(H, subMat, k, q - 1);
            
            r = min(k + 4, n);
            subMat = mat_mul_mat(mat_slice_mat(H, 0, r - 1, k, k + 2), mat_I_sub_beta_vvT(beta, v));
            mat_put_mat(H, subMat, 0, k);
            
            x = H[k + 1][k];
            y = H[k + 2][k];
            if (k < n - 3) {
                z = H[k + 3][k];
            }

            subMat = eye(n);
            mat_put_mat(subMat, mat_I_sub_beta_vvT(beta, v), k, k);
            
            P = mat_mul_mat(P, subMat);

        }
        u.resize(2); u[0] = x; u[1] = y;
        v.resize(2);
        beta = householder(u, v);
        
        subMat = mat_mul_mat(mat_I_sub_beta_vvT(beta, v), mat_slice_mat(H, n - 2, n - 1, n - 3, n - 1));
        mat_put_mat(H, subMat, n - 2, n - 3);
        
        subMat = mat_mul_mat(mat_slice_mat(H, 0, n - 1, n - 2, n - 1), mat_I_sub_beta_vvT(beta, v));
        mat_put_mat(H, subMat, 0, n - 2);

        subMat = eye(n);
        mat_put_mat(subMat, mat_I_sub_beta_vvT(beta, v), n - 2, n - 2);
        
        P = mat_mul_mat(P, subMat);

        for (int i = 2; i < n; i++) {
            for (int j = 0; j < i - 1; j++) {
                H[i][j] = 0;
            }
        }
    }
}

vector<vector<double> > transpose(vector<vector<double> > A) {
    // 返回矩阵A的转置
    int m = A.size();
    int n = A[0].size();
    vector<vector<double> > AT(n, vector<double>(m));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            AT[j][i] = A[i][j];
        }
    }
    return AT;
}

bool isApproximable(vector<vector<double> > H) {
    // 判断Hessenberg矩阵H是否可约
    int n = H.size();
    for (int i = 0; i < n - 1; i++) {
        if (H[i + 1][i] == 0) {
            return true;
        }
    }
    return false;
}

bool isProposedUpTriMat(vector<vector<double> > H) {
    // 判断Hessenberg矩阵H是否是拟上三角阵
    int n = H.size();
    double delta;
    if (n <= 1) {
        return true;
    }
    if (n == 2) {
        if (H[1][1] == 0) {
            return true;
        }
        else {
            if (pow(H[0][0] + H[1][1], 2) - 4 * (H[0][0] * H[1][1] - H[0][1] * H[1][0]) < 0) {
                return true;
            }
        }
        return false;
    }
    for (int i = 0; i < n - 1; i++) {
        if (H[i + 1][i]) {  // 当前考虑的次对角元非零
            if (i == 0) {
                if (H[i + 2][i + 1]) {  // 与之相邻的次对角元也非零
                    return false;
                }
            }
            else if (i == n - 2) {
                if (H[i][i - 1]) {  // 与之相邻的次对角元也非零
                    return false;
                }
            }
            else {
                if (H[i + 2][i + 1] || H[i][i - 1]) {   // 与之相邻的次对角元也非零
                    return false;
                }
            }
            delta = pow(H[i][i] + H[i + 1][i + 1], 2) 
                - 4 * (H[i][i] * H[i + 1][i + 1] - H[i][i + 1] * H[i + 1][i]);
            if (delta > 0) {
                return false;
            }
        }
    }
    return true;
}

void find_m_l(vector<vector<double> > H, int& m, int& l) {
    // 算法6.4.3的辅助函数，完成(3)(ii)的功能，即找到满足要求的最大的非负整数m和最小的非负整数l
    int n = H.size();
    vector<vector<double> > H33;
    for (m = n; m > 0; m--) {
        H33 = mat_slice_mat(H, n - m, n - 1, n - m, n - 1);
        if (m == n) {
            if (isProposedUpTriMat(H33)) {
                break;
            }
        }
        else {
            if (isProposedUpTriMat(H33)&&(H[n - m][n - m - 1] == 0)) {
                break;
            }
        }
    }
    if (m == n) {
        l = 0;
    }
    else {
        for (l = 0; l < n - m; l++) {
            if (!isApproximable(mat_slice_mat(H, l, n - m - 1, l, n - m - 1))) {
                break;
            }
        }
    }
    
}

void implicitQR(vector<vector<double> >& A, double u) {
    // 计算实矩阵A的实Schur分解：隐式QR算法，对应于教材算法6.4.3
    // 计算A的实Schur标准型，并保存在A中；u-机器精度
    int n = A.size();
    int m = 0, l = 0;
    vector<vector<double> > H22, H12, H23, P;
    upperHessenberg(A);
    
    for (int i = 1; i < n; i++) {
        if (fabs(A[i][i - 1]) < (fabs(A[i][i]) + fabs(A[i - 1][i - 1])) * u) {
            A[i][i - 1] = 0;
        }
    }
    
    find_m_l(A, m, l);

    int count = 0;

    while (m != n) {

        count++;

        if (n - l - m > 2) {
            H22 = mat_slice_mat(A, l, n - m - 1, l, n - m - 1);
            doubleShiftQR(H22, P);
            mat_put_mat(A, H22, l, l);
        }
        else {
            P = eye(2);
            double t = ((H22[0][0] - H22[1][1]) 
                + sqrt(pow(H22[0][0] - H22[1][1], 2) + 4 * H22[0][1] * H22[1][0])) / (2 * H22[0][1]);
            P[0][0] = 1.0 / sqrt(1 + pow(t, 2));
            P[0][1] = t / sqrt(1 + pow(t, 2));
            P[1][0] = -P[0][1];
            P[1][1] = P[0][0];
        }
        
        if (l > 0) {
            H12 = mat_slice_mat(A, 0, l - 1, l, n - m - 1);//l=0咋办
            H12 = mat_mul_mat(H12, P);
            mat_put_mat(A, H12, 0, l);
        }

        if (m > 0) {
            H23 = mat_slice_mat(A, l, n - m - 1, n - m, n - 1);
            H23 = mat_mul_mat(transpose(P), H23);
            mat_put_mat(A, H23, l, n - m);//m=0咋办
        }
        
        for (int i = 1; i < n; i++) {
            if (fabs(A[i][i - 1]) < (fabs(A[i][i]) + fabs(A[i - 1][i - 1])) * u) {
                A[i][i - 1] = 0;
            }
        }

        if (count > 5000) {
            cout << "Max iteration 5000 reached!" << endl;
            break;
        }

        find_m_l(A, m, l);
    }
    cout << "Implicit QR Iterations: " << count << endl;
}

void implicitQREigenvalue(vector<vector<double> > A, vector<double>& Re, vector<double>& Im) {
    // 使用隐式QR迭代，计算矩阵A的全部特征值
    // 计算出的特征值实部保存在向量Re中，虚部保存在向量Im中，按照下标一一对应
    int n = A.size();
    if (A[0].size() != n) {
        cout << "ERROR: The input matrix is not a square!" << endl;
        return;
    }
    double u = 1e-7;
    implicitQR(A, u);
    Re.clear();
    Im.clear();
    double b, c, delta;
    int i = 0;
    while (i < n) {
        if ((i == n - 1) || (A[i + 1][i] == 0)) {
            Re.push_back(A[i][i]);
            Im.push_back(0);
            i++;
        }
        else {
            b = -(A[i][i] + A[i + 1][i + 1]);
            c = A[i][i] * A[i + 1][i + 1] - A[i][i + 1] * A[i + 1][i];
            delta = pow(b, 2) - 4 * c;
            if (delta >= 0) {
                cout << "ERROR: real eigenvalue mistaked!" << endl;
                return;
            }
            Re.push_back(-b / 2);
            Im.push_back(sqrt(-delta) / 2);
            Re.push_back(-b / 2);
            Im.push_back(-sqrt(-delta) / 2);
            i += 2;
        }
    }
}



// =============================第七章中定义的函数=============================

double nonDiagNorm(vector<vector<double> > A) {
    // 计算方阵A的非对角“范数”，其定义见教材P211
    int n = A.size();
    double norm = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            norm += pow(A[i][j], 2);
        }
    }
    for (int i = 0; i < n; i++) {
        norm -= pow(A[i][i], 2);
    }
    norm = sqrt(norm);
    return norm;
}

bool isPassed(vector<vector<double> > A, double delta, int& p, int& q) {
    // 判断过关Jacobi方法中，当前矩阵是否“过关”
    // delta-当前的关值；p, q用于存储绝对值超过关值的非对角元所在平面，若不存在则置为-1
    int n = A.size();
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(A[i][j]) > delta) {
                p = i;
                q = j;
                return false;
            }
        }
    }
    p = -1;
    q = -1;
    return true;
}

void rotateJacobi(vector<vector<double> > A, int p, int q, double& c, double& s) {
    // 计算J(p,q,\theta)中的\cos\theta和\sin\theta，保存在c, s中
    int n = A.size();
    if (fabs(A[p][q]) < 1e-9) {
        c = 1;
        s = 0;
    }
    else {
        double tau = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        double t = -1;
        if (tau >= 0) {
            t = 1;
        }
        t /= fabs(tau) + sqrt(1 + pow(tau, 2));
        c = 1.0 / sqrt(1 + pow(t, 2));
        s = t * c;
    }
}

vector<vector<double> > leftJacobiMul(vector<vector<double> > A, int p, int q, double c, double s) {
    // 用Jacobi矩阵J(p,q,\theta)左乘矩阵A，返回计算结果
    vector<vector<double> > result = A;
    int n = A.size();
    for (int i = 0; i < n; i++) {
        result[p][i] = c * A[p][i] + s * A[q][i];
        result[q][i] = -s * A[p][i] + c * A[q][i];
    }
    return result;
}

vector<vector<double> > rightJacobiMul(vector<vector<double> > A, int p, int q, double c, double s) {
    // 用Jacobi矩阵J(p,q,\theta)右乘矩阵A，返回计算结果
    vector<vector<double> > result = A;
    int n = A.size();
    for (int i = 0; i < n; i++) {
        result[i][p] = c * A[i][p] - s * A[i][q];
        result[i][q] = s * A[i][p] + c * A[i][q];
    }
    return result;
}

void passingJacobiMethod(vector<vector<double> > A, vector<double>& EigenValues, vector<vector<double> >& EigenVectors, double sigma) {
    // 过关Jacobi方法求对称矩阵的全部特征值和特征向量
    // A-实对称矩阵；EigenValues-全部特征值；EigenVextors-对应的规范特征向量组成的矩阵
    double delta = nonDiagNorm(A);
    int n = A.size();
    int p, q, count = 0;
    vector<vector<double> > Q, tempJ;
    EigenValues.resize(n);
    Q = eye(n);
    double c, s;
    while (delta > 1e-7) {
        while (!isPassed(A, delta, p, q)) {
            rotateJacobi(A, p, q, c, s);
            A = leftJacobiMul(A, p, q, c, -s);
            A = rightJacobiMul(A, p, q, c, s);
            Q = rightJacobiMul(Q, p, q, c, s);
            count++;
        }
        delta /= sigma;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(Q[i][j]) < 1e-7) Q[i][j] = 0;
            if (i == j) {
                EigenValues[i] = A[i][i];
                continue;
            }
            A[i][j] = 0;
        }
    }
    EigenVectors = Q;
    // cout << "Iterations: " << count;
}

void print_mat(vector<vector<double> > A, int digits) {
    // 打印矩阵A，其中digits为小数点后位数（重载）
    int m = A.size();
    int n = A[0].size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << setiosflags(ios::fixed);
            cout << setprecision(digits) << A[i][j] << '\t';
        }
        cout << endl;
    }
}

void print_vec(vector<double> v, int digits) {
    // 打印向量v，其中digits为小数点后位数（重载）
    int n = v.size();
    for (int i = 0; i < n; i++) {
        cout << setiosflags(ios::fixed);
        cout << setprecision(digits) << v[i] << '\t';
    }
    cout << endl;
}

int changeSignNum(vector<vector<double> > A, double mu) {
    // 计算 对称三对角矩阵 A在mu处的变号数s_n(mu)
    int n = A.size();
    vector<double> x(n), y(n);
    x[0] = A[0][0];
    for (int i = 1; i < n; i++) {
        x[i] = A[i][i];
        y[i] = A[i - 1][i];
    }
    int s = 0;
    double q = x[0] - mu;
    for (int k = 0; k < n; k++) {
        if (q < 0) {
            s++;
        }
        if (k < n - 1) {
            if (q == 0) {
                q = fabs(y[k + 1]) * 1e-9;
            }
            q = x[k + 1] - mu - y[k + 1] * y[k + 1] / q;
        }
    }
    return s;
}

double bisectionMethod(vector<vector<double> > A, int m, double eps) {
    // 使用二分法，求矩阵A的第m小的特征值
    // 其中，矩阵A要求是实对称三对角矩阵，eps是计算精度
    int n = A.size();
    double u = mat_infty_norm(A);
    double l = -u;
    double s, r;
    int count = 0;
    while (u - l > eps) {
        r = (l + u) / 2;
        s = changeSignNum(A, r);
        if (s >= m) {
            u = r;
        }
        else {
            l = r;
        }
        count++;
    }
    cout << "Bisection Iterations: " << count << endl;
    return (l + u) / 2;
}

vector<double> inversePowerMethod(vector<vector<double> > A, double lambda, double eps) {
    // 使用反幂法求解矩阵A的特征值lambda对应的特征向量
    // eps为计算精度
    int n = A.size();
    vector<double> z(n);
    vector<int> u;
    double norm = 0;
    double pre = 1, cur = 0;
    z[0] = 1;
    for (int i = 0; i < n; i++) {
        A[i][i] -= lambda;
    }
    gauss_elim_col_pivoting(A, u);
    while (fabs(pre - cur) > eps) {
        pre = vec_2_norm(z);
        vector_pb(u, z);
        forward_subs1(A, z);
        back_subs(A, z);
        norm = vec_2_norm(z);
        for (int k = 0; k < n; k++) {
            z[k] /= norm;
        }
        cur = vec_2_norm(z);
    }
    return z;
}

void helper(vector<double> eigen) {
    // 打表格小助手
    int n = eigen.size();
    sort(eigen.begin(), eigen.end());
    for (int i = 0; i < n; i++) {
        cout << setiosflags(ios::fixed);
        cout << setprecision(4) << eigen[i];
        if (i % 10 == 9) cout << " \\\\" << endl;
        else cout << " & ";
    }
}



// =============================第八章中定义的函数=============================

vector<vector<double> > xyT(vector<double> x, vector<double> y) {
    // 计算矩阵 xy^T，其中x, y为行向量，大小分别为m, n
    int m = x.size(), n = y.size();
    vector<vector<double> > result(m, vector<double>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = x[i] * y[j];
        }
    }
    return result;
}

vector<double> numMul(double lambda, vector<double> v) {
    // 计算向量数乘，返回计算结果
    vector<double> result = v;
    for (int i = 0; i < v.size(); i++) {
        result[i] *= lambda;
    }
    return result;
}

vector<vector<double> > numMul(double lambda, vector<vector<double> > A) {
    // 计算矩阵数乘，返回计算结果
    vector<vector<double> > result = A;
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            result[i][j] *= lambda;
        }
    }
    return result;
}

vector<vector<double> > biDiagonalization(vector<vector<double> > A, vector<vector<double> >& U, vector<vector<double> >& V) {
    // 对矩阵 A 用 Householder 变换做二对角化
    // 对应教材P233算法7.6.1
    // U^T A V = [B; 0]
    int n = A[0].size(), m = A.size();
    double beta;
    vector<double> v, x;
    vector<vector<double> > P, B;
    U = eye(m);
    V = eye(n);
    for (int k = 0; k < n; k++) {
        x = mat_slice_vec(A, k, m - 1, k);
        v.resize(x.size());
        beta = householder(x, v);
        mat_put_mat(A, mat_mul_mat(mat_I_sub_beta_vvT(beta, v), mat_slice_mat(A, k, m - 1, k, n - 1)), k, k);
        P = eye(m);
        mat_put_mat(P, mat_I_sub_beta_vvT(beta, v), k, k);
        U = mat_mul_mat(P, U);
        if (k < n - 2) {
            x = mat_slice_vec(transpose(A), k + 1, n - 1, k);
            v.resize(x.size());
            beta = householder(x, v);
            mat_put_mat(A, mat_mul_mat(mat_slice_mat(A, k, m - 1, k + 1, n - 1), mat_I_sub_beta_vvT(beta, v)), k, k + 1);
            P = eye(n);
            mat_put_mat(P, mat_I_sub_beta_vvT(beta, v), k + 1, k + 1);
            V = mat_mul_mat(V, P);
        }
    }
    U = transpose(U);
    B = eye(n);
    for (int i = 0; i < n - 1; i++) {
        B[i][i] = A[i][i];
        B[i][i + 1] = A[i][i + 1];
    }
    B[n - 1][n - 1] = A[n - 1][n - 1];
    return B;
}

void givens(double& c, double& s, double a, double b) {
    // 找Givens变换中的\cos\theta和\sin\theta，分别保存在c, s中
    double t;
    if (b == 0) {
        c = 1;
        s = 0;
    }
    else {
        if (fabs(b) > fabs(a)) {
            t = a / b;
            s = 1.0 / sqrt(1 + pow(t, 2));
            c = s * t;
        }
        else {
            t = b / a;
            c = 1.0 / sqrt(1 + pow(t, 2));
            s = c * t;
        }
    }
}

vector<vector<double> > WilkinsonSVD(vector<vector<double> > B, vector<vector<double> >& P, vector<vector<double> >& Q) {
    // 带Wilkinson位移的SVD迭代
    // 对应于算法7.6.2，这里要求 n>=3，P, Q要求大小和B相符合
    // newB = P^T B Q
    int n = B.size();
    if (n < 3) {
        cout << "ERROR: n < 3!" << endl;
    }
    double alpha, beta, delta, mu, y, z, c, s;
    alpha = pow(B[n - 1][n - 1], 2) + pow(B[n - 2][n - 1], 2);
    delta = (pow(B[n - 2][n - 2], 2) + pow(B[n - 3][n - 2], 2) - alpha) / 2;
    beta = B[n - 2][n - 2] * B[n - 2][n - 1];
    mu = alpha - pow(beta, 2) / (delta + ((delta > 0) - (delta < 0)) * sqrt(pow(delta, 2) + pow(beta, 2)));
    y = pow(B[0][0], 2) - mu;
    z = B[0][0] * B[0][1];
    for (int k = 0; k < n - 1; k++) {
        givens(c, s, y, z);
        B = rightJacobiMul(B, k, k + 1, c, -s);
        Q = rightJacobiMul(Q, k, k + 1, c, -s);
        y = B[k][k];
        z = B[k + 1][k];
        givens(c, s, y, z);
        P = leftJacobiMul(P, k, k + 1, c, s);
        if (k < n - 2) {
            B = leftJacobiMul(B, k, k + 1, c, s);
            y = B[k][k + 1];
            z = B[k][k + 2];
        }
        else {
            B = leftJacobiMul(B, k, k + 1, c, s);
        }
    }
    P = transpose(P);
    return B;
}

void find_p_q(vector<vector<double> > B, int& p, int& q) {
    // 算法7.6.3中的辅助程序，找满足要求的p和q
    int n = B.size(), i;
    p = 0;
    q = 0;
    for (i = 0; i < n - 1; i++) {
        if (B[n - 2 - i][n - 1 - i] != 0) break;
    }
    if (i == n - 1) q = n;
    else {
        q = i;
        for (i = 1; i < n - q - 1; i++) {
            if (B[n - q - 2 - i][n - q - 1 - i] == 0) break;
        }
        p = n - q - 1 - i;
    }
}

void SVD(vector<vector<double> > A, vector<vector<double> >& U, vector<vector<double> >& Sigma, vector<vector<double> >& V, double eps) {
    // 矩阵A的SVD分解，其中A的大小为m*n，m\geq n
    // A = U^T \Sigma V
    int m = A.size();
    int n = A[0].size();
    int p, q, r, s;
    int count = 0;
    double cos, sin, a, b, c, d;
    if (m < n) {
        cout << "ERROR: m < n!" << endl;
        return;
    }
    vector<vector<double> > B, B22, H, IU, IV;
    U = eye(m); V = eye(n); // 初始化U, V
    B = biDiagonalization(A, U, V);

    U = transpose(U);

    double infNorm = mat_infty_norm(B);
    for (int i = 0; i < n - 1; i++) {
        if (fabs(B[i][i + 1]) <= (fabs(B[i][i]) + fabs(B[i + 1][i + 1])) * eps) {
            B[i][i + 1] = 0;
        }
    }
    for (int i = 0; i < n; i++) {
        if (fabs(B[i][i]) < infNorm * eps) {
            B[i][i] = 0;
        }
    }
    find_p_q(B, p, q);
    
    while (q < n) {
        s = 0;
        r = n - p - q;
        B22 = mat_slice_mat(B, p, p + r - 1, p, p + r - 1);
        for (int i = 0; i < r - 1; i++) {
            if (B22[i][i] == 0) {
                s++;
                for (int j = i + 1; j < r - 1; j++) {
                    givens(cos, sin, B22[j][j], B22[i][j]);
                    a = B[i][j]; c = B[j][j]; d = B[j][j + 1];
                    B[i][j] = -sin * a + cos * c;
                    B[i][j + 1] = cos * d;
                    B[j][j] = cos * a + sin * c;
                    B[j][j + 1] = sin * d;
                    H = U;
                    for (int k = 0; k < m; k++) {
                        H[p + i][k] = -sin * U[p + i][k] + cos * U[p + j][k];
                        H[p + j][k] = cos * U[p + i][k] + sin * U[p + j][k];
                    }
                    U = H;
                    H.clear();
                }
                int j = r - 1;
                givens(cos, sin, B22[j][j], B22[i][j]);
                a = B22[i][j]; b = B22[j][j];
                B22[j][j] = cos * a + sin * b;
                B22[i][j] = -sin * a + cos * b;
                H = U;
                for (int k = 0; k < m; k++) {
                    H[p + i][k] = -sin * U[p + i][k] + cos * U[p + j][k];
                    H[p + j][k] = cos * U[p + i][k] + sin * U[p + j][k];
                }
                U = H;
                H.clear();
            }
        }
        if (s == 0) {
            count++;
            if (r == 2) {
                a = B22[0][0]; b = B22[0][1]; d = B22[1][1];
                cos = (a + d) / sqrt(pow(a + d, 2) + pow(b, 2));
                sin = -b / sqrt(pow(a + d, 2) + pow(b, 2));
                B22 = leftJacobiMul(B22, 0, 1, cos, sin);
                U = leftJacobiMul(U, p, p + 1, cos, sin);
                vector<double> values;
                vector<vector<double> > vectors;
                passingJacobiMethod(B22, values, vectors, 10);
                B[p][p] = values[0]; B[p][p + 1] = 0;
                B[p + 1][p] = 0; B[p + 1][p + 1] = values[1];
                a = vectors[0][0]; b = vectors[0][1];
                c = vectors[1][0]; d = vectors[1][1];
                mat_put_mat(U, mat_mul_mat(transpose(vectors), mat_slice_mat(U, p, p + 1, 0, m - 1)), p, 0);
                mat_put_mat(V, mat_mul_mat(mat_slice_mat(V, 0, n - 1, p, p + 1), vectors), 0, p);
                a = (B22[0][0] > 0) - (B22[0][0] < 0);
                b = (B22[1][1] > 0) - (B22[1][1] < 0);
                H = { {a, 0}, {0, b} };
                mat_put_mat(U, mat_mul_mat(H, mat_slice_mat(U, p, p + 1, 0, m - 1)), p, 0);
                H.clear();
            }
            else {
                count++;
                IU = eye(r); IV = eye(r);
                B22 = WilkinsonSVD(B22, IU, IV);
                IU = transpose(IU);

                mat_put_mat(U, mat_mul_mat(IU, mat_slice_mat(U, p, p + r - 1, 0, m - 1)), p, 0);
                mat_put_mat(V, mat_mul_mat(mat_slice_mat(V, 0, n - 1, p, p + r - 1), IV), 0, p);
                mat_put_mat(B, B22, p, p);
            }
        }

        infNorm = mat_infty_norm(B);
        for (int i = 0; i < n - 1; i++) {
            if (fabs(B[i][i + 1]) <= (fabs(B[i][i]) + fabs(B[i + 1][i + 1])) * eps) {
                B[i][i + 1] = 0;
            }
        }
        for (int i = 0; i < n; i++) {
            if (fabs(B[i][i]) < infNorm * eps) {
                B[i][i] = 0;
            }
        }

        find_p_q(B, p, q);

    }

    for (int i = 0; i < n; i++) {
        if (B[i][i] < 0) {
            B[i][i] *= -1;
            mat_put_mat(V, mat_mul_mat(mat_slice_mat(V, 0, n - 1, i, i), { {-1} }), 0, i);
        }
    }

    V = transpose(V);

    vector<vector<double> > result(m, vector<double>(n));
    Sigma = result;
    for (int i = 0; i < n; i++) {
        Sigma[i][i] = B[i][i];
    }
    
    cout << "Iterations: " << count << endl;

}

double maxAbs(vector<vector<double> > A) {
    // 返回矩阵A中模最大的元素
    double max = -1;
    int m = A.size();
    int n = A[0].size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(A[i][j]) > max) {
                max = fabs(A[i][j]);
            }
        }
    }
    return max;
}