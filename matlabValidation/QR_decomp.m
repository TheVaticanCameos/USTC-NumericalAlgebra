A = [1 2;3 4;5 6];
m = 3;
n = 2;

for j=1:n
    if j<m
        [v, beta] = house(A(j:m,j));
        A(j:m,j:n) = (eye(m-j+1) - beta*(v')*v)*A(j:m,j:n);
        d(j) = beta;
        A(j+1:m,j) = v(2:m-j+1);
    end
end

