#include "Exercise.h"
clock_t start, ende;


// ==============================第一章习题==============================

void exercise_1_1()
{
	// initializations
	int N = 50;
	vector<vector<double> > A(N, vector<double>(N));
	vector<double> b(N);

	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = 6;
		A[i + 1][i] = 8;
		A[i][i + 1] = 1;
		b[i] = 15;
	}
	A[N - 1][N - 1] = 6;
	b[0] = 7;
	b[N - 1] = 14;

	int mode;
	cout << "Please choose the solving method:\n";
	cout << "1-Gaussian elimination\t 2-full pivoting Gaussian elimination\t 3-column pivoting Gaussian elimination" << endl;
	cin >> mode;

	// choose the solving method
	switch (mode) {
	case 1: {	// Gaussian elimination
		start = clock();
		gauss_elim(A);
		forward_subs1(A, b);
		back_subs(A, b);
		ende = clock();
		break; }
	case 2: {	// full pivoting Gaussian elimination
		vector<int> u, v;
		start = clock();
		gauss_elim_full_pivoting(A, u, v);
		vector_pb(u, b);
		forward_subs1(A, b);
		back_subs(A, b);
		vector_qb(v, b);
		ende = clock();
		break; }
	case 3: {	// column pivoting Gaussian elimination
		vector<int> u;
		start = clock();
		gauss_elim_col_pivoting(A, u);
		vector_pb(u, b);
		forward_subs1(A, b);
		back_subs(A, b);
		ende = clock();
		break; }
	default:
		break;
	}

	// print the solution
	cout << endl << "The solution is: " << endl;
	for (int i = 0; i < N; i++)
		cout << b[i] << " ";
	cout << endl << endl;

	// calculate and print the error
	double error = 0;
	for (int i = 0; i < N; i++) {
		error += (b[i] - 1) * (b[i] - 1);
	}
	error = sqrt(error);
	cout << "Error : " << error << endl;

	cout << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << "seconds";
}

void exercise_1_2_1()
{
	// initializations
	int N = 100;
	vector<vector<double> > A(N, vector<double>(N));
	vector<double> b(N);
	vector<double> x(N);
	for (int i = 0; i < N - 1; i++) {
		A[i][i] = 10;
		A[i][i + 1] = 1;
		A[i + 1][i] = 1;
		x[i] = rand() % 20;
	}
	A[N - 1][N - 1] = 10;
	x[N - 1] = rand() % 20;
	mat_mul_vec(A, x, b);

	// solve the equations using Cholesky decomposition
	int mode;
	cout << "Please choose the solving method:\n";
	cout << "1-Cholesky decomposition\t 2-modified Cholesky decomposition\t" << endl;
	cin >> mode;

	switch (mode){
	case 1: {
		start = clock();
		cholesky_decomp(A);
		forward_subs(A, b);
		back_subs(A, b);
		ende = clock();
		break;
	}
	case 2:{
		start = clock();
		modified_cholesky_decomp(A);
		forward_subs1(A, b);
		matrix_DLT(A);
		back_subs(A, b);
		ende = clock();
		break;
	}
	default:
		break;
	}
	
	double error = 0;
	for (int i = 0; i < N; i++) {
		error += (b[i] - x[i]) * (b[i] - x[i]);
	}
	error = sqrt(error);
	
	// print the solution
	cout << endl << "The solution is: " << endl;
	for (int i = 0; i < N; i++) {
		cout << setiosflags(ios::fixed);
		cout << setprecision(4) << b[i] << '\t';
		if ((i + 1) % 10 == 0)
			cout << endl;
	}

	cout << endl << "Error : " << error << endl;
	cout << endl << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << "seconds" << endl;
}

void exercise_1_2_2()
{
	// initialization
	int N = 40;
	vector<vector<double> > A(N, vector<double>(N));
	vector<double> b(N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = 1.0 / (i + j + 1);
			b[i] += A[i][j];
		}
	}

	// solve the equations using Cholesky decomposition
	int mode;
	cout << "Please choose the solving method:\n";
	cout << "1-Cholesky decomposition\t 2-modified Cholesky decomposition\t" << endl;
	cin >> mode;

	start = clock();
	switch (mode) {
	case 1: {
		cholesky_decomp(A);
		forward_subs(A, b);
		back_subs(A, b);
		break;
	}
	case 2: {
		modified_cholesky_decomp(A);
		forward_subs1(A, b);
		matrix_DLT(A);
		back_subs(A, b);
		break;
	}
	default:
		break;
	}

	
	// print the solution
	cout << endl;
	for (int i = 0; i < N; i++) {
		cout << setiosflags(ios::fixed);
		cout << setprecision(4) << b[i] << '\t';
		if ((i + 1) % 5 == 0)
			cout << endl;
	}
	ende = clock();

	double error = 0;
	for (int i = 0; i < N; i++) {
		error += pow(b[i] - 1, 2);
	}
	error = sqrt(error);
	cout << "Error: " << error << endl;
	// print the time consumed
	cout << endl << "Time consumed: " << (double)(ende - start) << " ms" << endl;
}

void exercise_1_3_1()
{
	// initializations
	int N = 100;
	vector<vector<double> > A(N, vector<double>(N));
	vector<double> b(N);
	for (int i = 0; i < N - 1; i++) {
		A[i][i] = 10;
		A[i][i + 1] = 1;
		A[i + 1][i] = 1;
		b[i] = rand() % 20;
	}
	A[N - 1][N - 1] = 10;
	b[N - 1] = rand() % 20;

	int mode;
	cout << "Please choose the solving method:\n";
	cout << "1-Gaussian elimination\t 2-full pivoting Gaussian elimination\t 3-column pivoting Gaussian elimination" << endl;
	cin >> mode;

	// choose the solving method
	switch (mode) {
	case 1: {	// Gaussian elimination
		start = clock();
		gauss_elim(A);
		forward_subs1(A, b);
		back_subs(A, b);
		ende = clock();
		break; }
	case 2: {	// full pivoting Gaussian elimination
		vector<int> u, v;
		start = clock();
		gauss_elim_full_pivoting(A, u, v);
		vector_pb(u, b);
		forward_subs1(A, b);
		back_subs(A, b);
		vector_qb(v, b);
		ende = clock();
		break; }
	case 3: {	// column pivoting Gaussian elimination
		vector<int> u;
		start = clock();
		gauss_elim_col_pivoting(A, u);
		vector_pb(u, b);
		forward_subs1(A, b);
		back_subs(A, b);
		ende = clock();
		break; }
	default:
		break;
	}

	// print the solution
	cout << "The solution is: " << endl;
	for (int i = 0; i < N; i++) {
		cout << setiosflags(ios::fixed);
		cout << setprecision(4) << b[i] << '\t';
		if ((i + 1) % 10 == 0)
			cout << endl;
	}

	cout << endl;
	cout << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << "seconds";
}

void exercise_1_3_2()
{
	// initialization
	int N = 40;
	vector<vector<double> > A(N, vector<double>(N));
	vector<double> b(N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = 1.0 / (i + j + 1);
			b[i] += A[i][j];
		}
	}

	int mode;
	cout << "Please choose the solving method:\n";
	cout << "1-Gaussian elimination\t 2-full pivoting Gaussian elimination\t 3-column pivoting Gaussian elimination" << endl;
	cin >> mode;

	// choose the solving method
	switch (mode) {
	case 1: {	// Gaussian elimination
		start = clock();
		gauss_elim(A);
		forward_subs1(A, b);
		back_subs(A, b);
		ende = clock();
		break; }
	case 2: {	// full pivoting Gaussian elimination
		vector<int> u, v;
		start = clock();
		gauss_elim_full_pivoting(A, u, v);
		vector_pb(u, b);
		forward_subs1(A, b);
		back_subs(A, b);
		vector_qb(v, b);
		ende = clock();
		break; }
	case 3: {	// column pivoting Gaussian elimination
		vector<int> u;
		start = clock();
		gauss_elim_col_pivoting(A, u);
		vector_pb(u, b);
		forward_subs1(A, b);
		back_subs(A, b);
		ende = clock();
		break; }
	default:
		break;
	}

	// print the solution
	cout << "The solution is: " << endl;
	for (int i = 0; i < N; i++) {
		cout << setiosflags(ios::fixed);;
		cout << setprecision(4) << b[i] << '\t';
		if ((i + 1) % 5 == 0)
			cout << endl;
	}
	
	cout << endl;
	cout << "Time consumed: " << (double)(ende - start) << " ms";
}



// ==============================第二章习题==============================

void exercise_2_1_1() {
	// 第二章习题1（1）
	for(int N = 5; N < 21; N++){
		vector<vector<double> > A(N, vector<double>(N));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				A[i][j] = 1.0 / (i + j + 1);
			}
		}
		double keppa = kappa_infty(A);
		cout << "N = " << N << '\t' << "kappa: " << keppa << endl;
	}
}

void exercise_2_1_2() {
	// 第二章习题1（2）
	int idx;
	double x_infty, error_infty;
	cout << "n \t" << "real precision \t" << "estimate precision: \t" << endl;
	for(int n = 5; n < 31; n++){
		vector<vector<double> > A(n, vector<double>(n));
		vector<double> x(n);
		vector<double> b(n);
		vector<double> sol(n);

		// 初始化矩阵A, 向量x
		srand((unsigned)time(NULL));
		for (int i = 0; i < n; i++) {
			A[i][i] = 1;
			A[i][n - 1] = 1;
			x[i] = rand() / (double)RAND_MAX * 50;	// x中的元素为0~50之间的浮点数
			for (int j = 0; j < i; j++) {
				A[i][j] = -1;
			}
		}

		mat_mul_vec(A, x, b);	// 计算向量b=Ax
		sol_eqs_col_LU(A, b, sol);	// 用列主元Gauss消去法求解线性方程组Ax=b，计算结果保存在sol中
		
		x_infty = vec_infty_norm(x, idx);
		for (int i = 0; i < n; i++) {
			x[i] -= sol[i];
		}
		error_infty = vec_infty_norm(x, idx);	// 计算解的实际精度

		cout << n << '\t' << error_infty / x_infty << '\t' << precision(A, b, sol) << endl;

	}
}



// ==============================第三章习题==============================

void exercise_3_1_1() {
	// 第三章习题1（1）
	int problem;
	cout << "Please choose the problem to solve: " << endl;
	cout << "1-chapter 1, exer 1    2-chapter 1, exer 2(1)    3-chapter 1, exer 2(2)" << endl;
	cin >> problem;
	switch (problem){
		case 1: {
			int N = 50;
			vector<vector<double> > A(N, vector<double>(N));
			vector<double> b(N);
			for (int i = 0; i < N - 1; i++)
			{
				A[i][i] = 6;
				A[i + 1][i] = 8;
				A[i][i + 1] = 1;
				b[i] = 15;
			}
			A[N - 1][N - 1] = 6;
			b[0] = 7;
			b[N - 1] = 14;

			start = clock();
			QR_solver(A, b);
			ende = clock();

			cout << endl << "The solution is: " << endl;
			print_vec(b);

			double error = 0;
			for (int i = 0; i < N; i++) {
				error += (b[i] - 1) * (b[i] - 1);
			}
			error = sqrt(error);
			cout << endl << "Error : " << error << endl;
			cout << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << " seconds";
			break;
		}
		case 2: {
			int N = 100;
			vector<vector<double> > A(N, vector<double>(N));
			vector<double> x(N);
			for (int i = 0; i < N - 1; i++) {
				A[i][i] = 10;
				A[i][i + 1] = 1;
				A[i + 1][i] = 1;
				x[i] = rand() % 20;
			}
			A[N - 1][N - 1] = 10;
			x[N - 1] = rand() % 20;
			vector<double> b(N);
			mat_mul_vec(A, x, b);

			start = clock();
			QR_solver(A, b);
			ende = clock();

			double error = 0;
			for (int i = 0; i < N; i++) {
				error += (b[i] - x[i]) * (b[i] - x[i]);
			}
			error = sqrt(error);

			cout << endl << "The solution is: " << endl;
			print_vec(b);
			cout << endl << "Error : " << error << endl;
			cout << endl << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << " seconds";
			break;
		}
		case 3: {
			int N = 40;
			vector<vector<double> > A(N, vector<double>(N));
			vector<double> b(N);
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					A[i][j] = 1.0 / (i + j + 1);
					b[i] += A[i][j];
				}
			}
			
			start = clock();
			QR_solver(A, b);
			ende = clock();

			cout << endl << "The solution is: " << endl;
			print_vec(b);

			double error = 0;
			for (int i = 0; i < N; i++) {
				error += pow(b[i] - 1, 2);
			}
			error = sqrt(error);
			cout << endl << "Error : " << error << endl;
			cout << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << " seconds";
			break;
		}
		default:
			break;
	}
}

void exercise_3_1_2() {
	// 第三章习题1（2）
	vector<double> t = {-1, -0.75, -0.5, 0, 0.25, 0.5, 0.75};
	vector<double> y = {1, 0.8125, 0.75, 1, 1.3125, 1.75, 2.3125};
	vector<vector<double> > A(7, vector<double>(3));
	for (int i = 0; i < 7; i++) {
		A[i][0] = t[i] * t[i];
		A[i][1] = t[i];
		A[i][2] = 1;
	}
	vector<double> coef(3);

	start = clock();
	double min = QR_LS(A, y, coef);
	ende = clock();
	
	cout << "optimal coef: " << endl;
	cout << "a\tb\tc" << endl;
	print_vec(coef);
	cout << "2 norm of the error vector: " << min << endl;
	cout << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << " seconds";
}

void exercise_3_1_3() {
	// 第三章习题1（3）
	vector<vector<double>> A =
	{ {1,4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0},
	{1,5.0208, 1, 3.531, 1.5, 2, 7, 4, 62, 1, 1, 0},
	{1,4.5429, 1, 2.275, 1.175, 1, 6, 3, 40,  2, 1, 0},
	{1,4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1, 0},
	{1,5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0},
	{1,3.891, 1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0},
	{1,5.898, 1, 5.85, 1.24, 1, 7, 3, 51, 2, 1,  1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,15.4202, 2.5,  9.8, 3.42, 2, 10, 5, 42, 2, 1, 1},
	{1,14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1},
	{1,5.8282, 1, 6.435, 1.225, 2, 6, 3, 32, 1, 1, 0},
	{1,5.3003, 1, 4.9883, 1.552, 1, 6, 3, 30, 1, 2, 0},
	{1,6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2, 0},
	{1,5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0},
	{1,5.05, 1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4, 1, 0},
	{1,6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1},
	{1,7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0},
	{1,9.0384, 1, 7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0},
	{1,5.9894, 1, 5.52, 1.256, 2, 6, 3, 40, 4, 1, 1},
	{1,7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1, 0},
	{1,8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1},
	{1,6.0931, 1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0},
	{1,8.3607, 1.5, 9.15, 1.777, 2., 8, 4, 48, 1, 1, 1},
	{1,8.14, 1, 8, 1.504, 2, 7, 3, 3, 1, 3, 0},
	{1,9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1, 0},
	{1,12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
	vector<double> b =
	{ 25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
	28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
	30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
	37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0 };
	vector<double> x(12);
	
	start = clock();
	double min = QR_LS(A, b, x);
	ende = clock();

	cout << "optimal coef: " << endl;
	for (int i = 0; i < 12; i++) {
		cout << "x" << i << '\t';
	}
	cout << endl;
	for (int i = 0; i < 12; i++) {
		cout << setiosflags(ios::fixed);
		cout << setprecision(4) << x[i] << '\t';
	}
	cout << endl << "2 norm of the error vector: " << min << endl;
	cout << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << " seconds";
}



// ==============================第四章习题==============================

void exercise_4_1() {
	// 第四章习题1
	int n = 100;
	double eps = 1.0;
	double a = 0.5;
	double h = 1.0 / n;
	double omega = 1.7;
	vector<vector<double> > A(n - 1, vector<double>(n - 1));
	vector<double> b(n - 1, a * h * h);	// 先用于存储线性方程组的右端项，求解完成后用于存储精确解
	vector<double> yCur(n - 1), yPre(n - 1, 1);	// yCur是迭代1步后的解，yPre是迭代1步前的解，迭代的初始值为1
	for (int i = 0; i < n - 2; i++) {
		A[i][i] = -(2 * eps + h);
		A[i][i + 1] = eps + h;
		A[i + 1][i] = eps;
	}
	A[n - 2][n - 2] = -(2 * eps + h);
	b[n - 2] = b[n - 2] - eps - h;

	int iterNum = 1;
	int option;
	cout << "epsilon: " << eps << endl;
	cout << "Please choose the solving method: " << endl;
	cout << "1-Jacobi iteration\n2-Gauss-Seidel iteration\n3-SOR iteration" << endl;
	cin >> option;
	switch (option){
		case 1:	// Jacobi迭代
			start = clock();
			Jacobi_preparation(A, b);
			yCur = Jacobi_iteration(yPre, A, b);
			while (vec_2_norm(vec_sub(yCur, yPre)) > 1e-6) {
				yPre = yCur;
				yCur = Jacobi_iteration(yPre, A, b);
				iterNum += 1;
			}
			ende = clock();
			break;
		case 2:	// Gauss-Seidel迭代
			start = clock();
			Jacobi_preparation(A, b);
			yCur = Gauss_Seidel_iteration(yPre, A, b);
			while (vec_2_norm(vec_sub(yCur, yPre)) > 1e-6) {
				yPre = yCur;
				yCur = Gauss_Seidel_iteration(yPre, A, b);
				iterNum += 1;
			}
			ende = clock();
			break;
		case 3:	// SOR迭代
			start = clock();
			yCur = SOR_iteration(yPre, A, b, omega);
			while (vec_2_norm(vec_sub(yCur, yPre)) > 1e-6) {
				yPre = yCur;
				yCur = SOR_iteration(yPre, A, b, omega);
				iterNum += 1;
			}
			ende = clock();
			break;
		default:
			break;
	}

	cout << "solution: " << endl;
	for (int i = 0; i < n - 1; i++) {
		cout << setiosflags(ios::fixed);
		cout << setprecision(4) << yCur[i] << '\t';
		if (i % 11 == 10) {
			cout << endl;
		}
	}

	cout << endl << "iteration times: " << iterNum << endl;
	cout << "time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << " s" << endl;

	double x;
	for (int i = 0; i < n - 1; i++) {
		x = (i + 1) * h;
		b[i] = (1 - a) * (1 - exp(-x / eps)) / (1 - exp(-1 / eps)) + a * x;
	}

	cout << "Error in 2 norm: " << vec_2_norm(vec_sub(yCur, b)) << endl;
	
}

void exercise_4_2() {
	// 第四章习题2
	int n = 80;
	double h = 1.0 / n;
	double eps = 1e-7;
	double omega = 1.92;
	int option;
	vector<vector<double> > sol;
	cout << "Please choose the solving method: " << endl;
	cout << "1-Jacobi iteration\n2-Gauss-Seidel iteration\n3-SOR iteration" << endl;
	cin >> option;
	switch (option){
		case 1:	// Jacobi迭代法
			start = clock();
			sol = pde2D_Jacobi(n, eps);
			ende = clock();
			break;
		case 2:	// Gauss-Seidel迭代法
			start = clock();
			sol = pde2D_GS(n, eps);
			ende = clock();
			break;
		case 3:	// SOR迭代法
			cout << "omega: " << omega << endl;
			start = clock();
			sol = pde2D_SOR(n, eps, omega);
			ende = clock();
			break;
		default:
			break;
	}

	int ii, jj;
	double min = find_min(sol, ii, jj);
	cout << "time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << " s" << endl;
	cout << "min: " << min << endl;
	cout << "x: " << (ii + 1) * h << endl;
	cout << "y: " << (jj + 1) * h << endl;
}



// ==============================第五章习题==============================

void exercise_5_1() {
	// 第五章习题1
	int n = 20;
	int m = pow(n - 1, 2);
	double h = 1.0 / n;
	vector<vector<double> > A(m, vector<double>(m));
	for (int i = 0; i < m; i++) {
		A[i][i] = 1 + pow(h, 2) / 4;
		if (i % (n - 1) != n - 2) {
			A[i][i + 1] = -0.25;
			A[i + 1][i] = -0.25;
		}
	}
	for (int i = 0; i < n - 2; i++) {
		for (int j = 0; j < n - 1; j++) {
			A[i * (n - 1) + j][(i + 1) * (n - 1) + j] = -0.25;
			A[(i + 1) * (n - 1) + j][i * (n - 1) + j] = -0.25;
		}
	}
	vector<double> b(m);
	for (int j = 0; j < n - 1; j++) {
		for (int i = 0; i < n - 1; i++) {
			b[j * (n - 1) + i] = pow(h, 2)/4 * sin((i + 1) * (j + 1) * h * h);
			if (i == 0) {
				b[j * (n - 1) + i] += 0.25 * (j + 1) * (j + 1) * h * h;
			}
			if (i == n - 2) {
				b[j * (n - 1) + i] += 0.25 * (1 + pow((j + 1) * h, 2));
			}
			if (j == 0) {
				b[j * (n - 1) + i] += 0.25 * (i + 1) * (i + 1) * h * h;
			}
			if (j == n - 2) {
				b[j * (n - 1) + i] += 0.25 * (1 + pow((i + 1) * h, 2));
			}
		}
	}
	vector<double> x(m);
	start = clock();
	CGD(A, x, b, 1e-7, m);
	ende = clock();
	cout << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << "seconds" << endl;
	cout << "solution: " << endl;
	for (int i = 0; i < m; i++) {
		cout << setiosflags(ios::fixed);
		cout << setprecision(4) << x[i] << '\t';
		if (i % 15 == 14) {
			cout << endl;
		}
	}	
}

void exercise_5_2() {
	// 第五章习题2
	int n = 60;
	vector<vector<double> > A(n, vector<double>(n));
	vector<double> b(n), x(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = 1.0 / (i + j + 1);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			b[i] += A[i][j];
		}
		b[i] /= 3;
	}
	start = clock();
	CGD(A, x, b, 1e-7, n);
	ende = clock();
	cout  << "Time consumed: " << (double)(ende - start) / CLOCKS_PER_SEC << "seconds" << endl;
	cout << "solution: " << endl;
	for (int i = 0; i < n; i++) {
		cout << setiosflags(ios::fixed);
		cout << setprecision(4) << x[i] << '\t';
		if ((i + 1) % 10 == 0) {
			cout << endl;
		}
	}
}

void exercise_5_3() {
	// 第五章习题3
	vector<vector<double> > A(5, vector<double>(5));
	vector<double> b(5), x(5), xCur(5), xPre(5);
	double eps = 1e-7;
	int iterNum = 0;
	A[0][0] = 10; A[0][1] = 1; A[0][2] = 2; A[0][3] = 3; A[0][4] = 4;
	A[1][0] = 1; A[1][1] = 9; A[1][2] = -1; A[1][3] = 2; A[1][4] = -3;
	A[2][0] = 2; A[2][1] = -1; A[2][2] = 7; A[2][3] = 3; A[2][4] = -5;
	A[3][0] = 3; A[3][1] = 2; A[3][2] = 3; A[3][3] = 12; A[3][4] = -1;
	A[4][0] = 4; A[4][1] = -3; A[4][2] = -5; A[4][3] = -1; A[4][4] = 15;
	b[0] = 12; b[1] = -27; b[2] = 14; b[3] = -17; b[4] = 12;

	int method;
	cout << "Please choose the solving method: " << endl;
	cout << "(1) Jacobi iteration\t (2) Gauss-Seidel iteration\t (3) Conjugate gradient" << endl;
	cin >> method;
	switch (method){
		case 1:
			start = clock();
			Jacobi_preparation(A, b);
			xCur = Jacobi_iteration(xPre, A, b);
			while (vec_2_norm(vec_sub(xCur, xPre)) > eps) {
				xPre = xCur;
				xCur = Jacobi_iteration(xPre, A, b);
				iterNum += 1;
			}
			ende = clock();
			x = xCur;
			break;
		case 2:
			start = clock();
			Jacobi_preparation(A, b);
			xCur = Gauss_Seidel_iteration(xPre, A, b);
			while (vec_2_norm(vec_sub(xCur, xPre)) > eps) {
				xPre = xCur;
				xCur = Gauss_Seidel_iteration(xPre, A, b);
				iterNum += 1;
			}
			ende = clock();
			x = xCur;
			break;
		case 3:
			start = clock();
			CGD(A, x, b, eps, 20);
			ende = clock();
			break;
		default:
			break;
	}

	if (method != 3) {
		cout << "iteration: " << iterNum << endl;
	}
	cout << "Time consumed: " << (double)(ende - start) << " ms" << endl;
	cout << "solution: " << endl;
	print_vec(x);
}



// ==============================第六章习题==============================

void exercise_6_1() {
	// 第六章习题1
	vector<double> coef;
	int num;
	cout << "Please choose the question number: ";
	cin >> num;
	switch (num){
		case 1:
			coef = { 3, -5, 1 };
			break;
		case 2:
			coef = { -1, -3, 0 };
			break;
		case 3:
			coef = { -1000, 790, -99902, 79108.9, 9802.08, 10891.01, 208.01, 101 };
			break;
		default:
			break;
	}
	double root;
	start = clock();
	root= maxAbsRoot(coef);
	ende = clock();
	cout << "Time consumed: " << (double)(ende - start) << " ms" << endl;
	cout << "Root with max abs: " << root << endl;
}

void exercise_6_2() {
	// 第六章习题2
	int num;
	cout << "Please choose the question number: ";
	cin >> num;
	double x;
	vector<vector<double> > coefMat(41, vector<double>(41));
	vector<vector<double> > A = { {9.1, 3.0, 2.6, 4.0}, 
								  {4.2, 5.3, 4.7, 1.6}, 
							      {3.2, 1.7, 9.4, 0.0}, 
								  {6.1, 4.9, 3.5, 6.2} };
	vector<double> Re, Im;
	for (int i = 0; i < 40; i++) {
		coefMat[i + 1][i] = 1;
	}
	coefMat[0][40] = -1; coefMat[0][37] = -1;
	switch (num) {
		case 1:
			start = clock();
			implicitQREigenvalue(coefMat, Re, Im);
			ende = clock();
			cout << "Time consumed: " << (double)(ende - start) << " ms" << endl;
			cout << "All roots: " << endl;
			for (int i = 0; i < 41; i++) {
				cout << setiosflags(ios::fixed);
				cout << setprecision(4) << Re[i] << " + " << Im[i] << " i";
				if (i % 4 == 3) {
					cout << endl;
				}
				else {
					cout << '\t';
				}
			}
			break;
		case 2:
			cout << "Please input the value of x: ";
			cin >> x;
			A[2][3] = x;
			start = clock();
			implicitQREigenvalue(A, Re, Im);
			ende = clock();
			cout << "Time consumed: " << (double)(ende - start) << " ms" << endl;
			cout << "All eigenvalues: " << endl;
			for (int i = 0; i < 4; i++) {
				cout << setiosflags(ios::fixed);
				cout << setprecision(4) << Re[i] << " + " << Im[i] << " i" << endl;
			}
			break;
		default:
			break;
	}
}



// ==============================第七章习题==============================

void exercise_7_1() {
	// 第七章习题1
	int n;
	double sigma = 10;
	cout << "Please input the size of matrix:" << endl;
	cin >> n;
	vector<vector<double> > A(n, vector<double>(n));
	for (int i = 0; i < n - 1; i++) {
		A[i][i] = 4;
		A[i][i + 1] = 1;
		A[i + 1][i] = 1;
	}
	A[n - 1][n - 1] = 4;
	
	cout << "n = " << n << ", ";

	vector<vector<double> > EigenVectors;
	vector<double> EigenValues;

	start = clock();
	passingJacobiMethod(A, EigenValues, EigenVectors, sigma);
	ende = clock();

	cout << ", " << "Time consumed: "<< (double)(ende - start) << " ms" << endl;
	cout << "Eigenvalues: " << endl;
	print_vec(EigenValues, 4);
	cout << endl;
	cout << "Eigenvectors(Q_k): " << endl;
	print_mat(EigenVectors, 4);

	/*cout << endl << endl;
	helper(EigenValues);*/
}

void exercise_7_2() {
	// 第七章习题2
	int n = 100;
	vector<vector<double> > A(n, vector<double>(n));
	for (int i = 0; i < n - 1; i++) {
		A[i][i] = 2;
		A[i][i + 1] = -1;
		A[i + 1][i] = -1;
	}
	A[n - 1][n - 1] = 2;
	
	start = clock();
	double maxEigenvalue = bisectionMethod(A, n, 1e-7);
	vector<double> maxV = inversePowerMethod(A, maxEigenvalue, 1e-7);
	ende = clock();
	cout << "Time consumed: " << (double)(ende - start) << " ms" << endl;
	cout << "Max Eigenvalue: " << maxEigenvalue << endl;
	cout << "Corresponding Eigenvector: " << endl;
	print_vec(maxV,  4);

	cout << endl;

	start = clock();
	double minEigenvalue = bisectionMethod(A, 1, 1e-7);
	vector<double> minV = inversePowerMethod(A, minEigenvalue, 1e-7);
	ende = clock();
	cout << "Time consumed: " << (double)(ende - start) << " ms" << endl;
	cout << "Min Eigenvalue: " << minEigenvalue << endl;
	cout << "Corresponding Eigenvector: " << endl;
	print_vec(minV, 4);
}



// ==============================第八章习题==============================

void exercise_8() {
	// 第八章习题
	vector<vector<double>> A = { {1.0000000000,4.9176000000,1.0000000000,3.4720000000,0.9980000000,1.0000000000,7.0000000000,4.0000000000,42.0000000000,3.0000000000,1.0000000000,0.0000000000},{1.0000000000,5.0208000000,1.0000000000,3.5310000000,1.5000000000,2.0000000000,7.0000000000,4.0000000000,62.0000000000,1.0000000000,1.0000000000,0.0000000000},{1.0000000000,4.5429000000,1.0000000000,2.2750000000,1.1750000000,1.0000000000,6.0000000000,3.0000000000,40.0000000000,2.0000000000,1.0000000000,0.0000000000},{1.0000000000,4.5573000000,1.0000000000,4.0500000000,1.2320000000,1.0000000000,6.0000000000,3.0000000000,54.0000000000,4.0000000000,1.0000000000,0.0000000000},{1.0000000000,5.0597000000,1.0000000000,4.4550000000,1.1210000000,1.0000000000,6.0000000000,3.0000000000,42.0000000000,3.0000000000,1.0000000000,0.0000000000},{1.0000000000,3.8910000000,1.0000000000,4.4550000000,0.9880000000,1.0000000000,6.0000000000,3.0000000000,56.0000000000,2.0000000000,1.0000000000,0.0000000000},{1.0000000000,5.8980000000,1.0000000000,5.8500000000,1.2400000000,1.0000000000,7.0000000000,3.0000000000,51.0000000000,2.0000000000,1.0000000000,1.0000000000},{1.0000000000,5.6039000000,1.0000000000,9.5200000000,1.5010000000,0.0000000000,6.0000000000,3.0000000000,32.0000000000,1.0000000000,1.0000000000,0.0000000000},{1.0000000000,15.4202000000,2.5000000000,9.8000000000,3.4200000000,2.0000000000,10.0000000000,5.0000000000,42.0000000000,2.0000000000,1.0000000000,1.0000000000},{1.0000000000,14.4598000000,2.5000000000,12.8000000000,3.0000000000,2.0000000000,9.0000000000,5.0000000000,14.0000000000,4.0000000000,1.0000000000,1.0000000000},{1.0000000000,5.8282000000,1.0000000000,6.4350000000,1.2250000000,2.0000000000,6.0000000000,3.0000000000,32.0000000000,1.0000000000,1.0000000000,0.0000000000},{1.0000000000,5.3003000000,1.0000000000,4.9883000000,1.5520000000,1.0000000000,6.0000000000,3.0000000000,30.0000000000,1.0000000000,2.0000000000,0.0000000000},{1.0000000000,6.2712000000,1.0000000000,5.5200000000,0.9750000000,1.0000000000,5.0000000000,2.0000000000,30.0000000000,1.0000000000,2.0000000000,0.0000000000},{1.0000000000,5.9592000000,1.0000000000,6.6660000000,1.1210000000,2.0000000000,6.0000000000,3.0000000000,32.0000000000,2.0000000000,1.0000000000,0.0000000000},{1.0000000000,5.0500000000,1.0000000000,5.0000000000,1.0200000000,0.0000000000,5.0000000000,2.0000000000,46.0000000000,4.0000000000,1.0000000000,1.0000000000},{1.0000000000,5.6039000000,1.0000000000,9.5200000000,1.5010000000,0.0000000000,6.0000000000,3.0000000000,32.0000000000,1.0000000000,1.0000000000,0.0000000000},{1.0000000000,8.2464000000,1.5000000000,5.1500000000,1.6640000000,2.0000000000,8.0000000000,4.0000000000,50.0000000000,4.0000000000,1.0000000000,0.0000000000},{1.0000000000,6.6969000000,1.5000000000,6.0920000000,1.4880000000,1.5000000000,7.0000000000,3.0000000000,22.0000000000,1.0000000000,1.0000000000,1.0000000000},{1.0000000000,7.7841000000,1.5000000000,7.1020000000,1.3760000000,1.0000000000,6.0000000000,3.0000000000,17.0000000000,2.0000000000,1.0000000000,0.0000000000},{1.0000000000,9.0384000000,1.0000000000,7.8000000000,1.5000000000,1.5000000000,7.0000000000,3.0000000000,23.0000000000,3.0000000000,3.0000000000,0.0000000000},{1.0000000000,5.9894000000,1.0000000000,5.5200000000,1.2560000000,2.0000000000,6.0000000000,3.0000000000,40.0000000000,4.0000000000,1.0000000000,1.0000000000},{1.0000000000,7.5422000000,1.5000000000,4.0000000000,1.6900000000,1.0000000000,6.0000000000,3.0000000000,22.0000000000,1.0000000000,1.0000000000,0.0000000000},{1.0000000000,8.7951000000,1.5000000000,9.8900000000,1.8200000000,2.0000000000,8.0000000000,4.0000000000,50.0000000000,1.0000000000,1.0000000000,1.0000000000},{1.0000000000,6.0931000000,1.5000000000,6.7265000000,1.6520000000,1.0000000000,6.0000000000,3.0000000000,44.0000000000,4.0000000000,1.0000000000,0.0000000000},{1.0000000000,8.3607000000,1.5000000000,9.1500000000,1.7770000000,2.0000000000,8.0000000000,4.0000000000,48.0000000000,1.0000000000,1.0000000000,1.0000000000},{1.0000000000,8.1400000000,1.0000000000,8.0000000000,1.5040000000,2.0000000000,7.0000000000,3.0000000000,3.0000000000,1.0000000000,3.0000000000,0.0000000000},{1.0000000000,9.1416000000,1.5000000000,7.3262000000,1.8310000000,1.5000000000,8.0000000000,4.0000000000,31.0000000000,4.0000000000,1.0000000000,0.0000000000},{1.0000000000,12.0000000000,1.5000000000,5.0000000000,1.2000000000,2.0000000000,6.0000000000,3.0000000000,30.0000000000,3.0000000000,1.0000000000,1.0000000000} };
	int m = A.size(), n = A[0].size();
	vector<vector<double> > U, Sigma, V;
	double eps = 1e-7;
	double eU, eV, eS;
	start = clock();
	SVD(A, U, Sigma, V, eps);	// A = U^T Sigma V
	ende = clock();
	cout << "Time consumed: " << (double)(ende - start) << " ms" << endl;
	vector<double> singular(n), copy;
	for (int i = 0; i < n; i++) {
		singular[i] = Sigma[i][i];
	}
	copy = singular;
	sort(copy.begin(), copy.end());
	cout << endl << "Singular values (increasing order): " << endl;
	print_vec(copy, 4);
	cout << endl;
	eU = maxAbs(mat_sub(mat_mul_mat(transpose(U), U), eye(m)));
	eV = maxAbs(mat_sub(mat_mul_mat(V, transpose(V)), eye(n)));
	eS = maxAbs(mat_sub(mat_mul_mat(transpose(U), mat_mul_mat(Sigma, V)), A));
	cout << setiosflags(ios::fixed);
	cout << "eU: " << setprecision(16) << eU << endl;
	cout << "eV: " << setprecision(16) << eV << endl;
	cout << "eS: " << setprecision(16) << eS << endl;
	cout << endl << "A = U^T Sigma V" << endl << endl;
	cout << "U: " << endl;
	print_mat(U, 4);
	cout << endl;
	cout << "V: " << endl;
	print_mat(V, 4);
	cout << endl;
	cout << "Sigma: " << endl;
	print_mat(Sigma, 4);

}