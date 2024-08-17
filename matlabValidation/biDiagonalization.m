function A = biDiagonalization(A)
    n = size(A, 2);
    m = size(A, 1);
    for k = 1:n
        [v, beta] = house(A(k:m, k));
        uT = beta*v*A(k:m, k:n);
        A(k:m, k:n) = A(k:m, k:n) - v'*uT;
        A(k+1:m, k) = v(2:m-k+1);
        b(k) = beta;
        if k < n-1
            [v, beta] = house(A(k, k+1:n)');
            u = beta*A(k:m, k+1:n)*v';
            A(k:m, k+1:n) = A(k:m, k+1:n) - u*v;
            A(k, k+2:n) = v(2:n-k)';
            c(k) = beta;
        end
    end
end