A = [1.1908, -1.0565, -2.1707, 0.5913, 0.0000, 0.7310;
    -1.2025, 1.4151, -0.0592, -0.6436, -0.3179, 0.5779;
    -0.0198, -0.8051, -1.0106, 0.3803, 1.0950, 0.0403;
    -0.1567, 0.5287, 0.6145, -1.0091, -1.8740, 0.6771;
    -1.6041, 0.2193, 0.5077, -0.0195, 0.4282, 0.5689;
    0.2573, -0.9219, 1.6924, -0.0482, 0.8956, -0.2556];
% A = [1, 0, 0, 0;
%     0, 1, -1, 0;
%     0, 1, 1, 0;
%     0, 0, 0, 1];
n = size(A, 1);
for k = 1: n-2
    [v, beta] = house(A(k+1:n, k));
    A(k+1:n, k:n) = (eye(length(v)) - beta*(v'*v))*A(k+1:n, k:n);
    A(1:n, k+1:n) = A(1:n, k+1:n)*(eye(length(v)) - beta*(v'*v));
end
% A(abs(A)<1e-7) = 0
% B = A;
% for i = 1:10
%     B = doubleShiftQR(B);
% end
% B(abs(B)<1e-7) = 0;
% for i=1:10
%     B = doubleShiftQR(B);
%     B(abs(B)<1e-7) = 0;
% end
