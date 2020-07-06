function [x,resvec] = mygmres(A,b,tol,maxits)
% MYGMRES 

if isnumeric(A)
    A = @(x) A*x;
end
resvec = norm(b);
V = b/resvec;
j = 1;
while j < maxits
    w = A(V(:,j));
    for i = 1:j
        H(i,j) = V(:,i)'*w;
        w = w - H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(w);
    V(:,j+1) = w/H(j+1,j);
    coeffs = resvec(1)*(H\eye(j+1,1));
    x = V(:,1:j)*coeffs;
    resvec(j+1) = norm(H*coeffs - eye(j+1,1)*resvec(1)); % = norm(b-A(x))
    
    if resvec(j+1) <= tol*resvec(1)
        return
    end
    j = j+1;
end