function H = doubleShiftQR(H)
    n = size(H, 1);
    m = n-1;
    s = H(m,m) + H(n,n);
    t = H(m,m)*H(n,n) - H(m,n)*H(n,m);
    x = H(1,1)*H(1,1) + H(1,2)*H(2,1) - s*H(1,1) + t;
    y = H(2,1)*(H(1,1) + H(2,2) - s);
    z = H(2,1)*H(3,2);
    for k = 0: n-3
        [v, beta] = house([x,y,z]');
        q = max(1,k);
        H(k+1:k+3, q:n) = (eye(length(v)) - beta*(v'*v))*H(k+1:k+3, q:n);
        r = min(k+4, n);
        H(1:r, k+1:k+3) = H(1:r, k+1:k+3)*(eye(length(v)) - beta*(v'*v));
        x = H(k+2, k+1);
        y = H(k+3, k+1);
        if k < n-3
            z = H(k+4, k+1);
        end
    end
    [v, beta] = house([x, y]');
    H(n-1:n, n-2:n) = (eye(length(v)) - beta*(v'*v))*H(n-1:n, n-2:n);
    H(1:n, n-1:n) = H(1:n, n-1:n)*(eye(length(v)) - beta*(v'*v));
end