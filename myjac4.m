function J =  myjac4(t,M1,M2,JFu,JFv,JGu,JGv,C,C2,u0,vT,uv)
%MYJAC4 
% builds the large Jacobian DR (faster version)
% this one uses the interlaced fwd form of the residual
% Note that C2 is not used.

dim = size(u0,1); n = length(t)-1;
UV = reshape(uv,2*dim,n+1);
u = UV(1:dim,:); v = UV(dim+1:end,:);

J = spalloc(2*dim*(n+1),2*dim*(n+1),round(0.1*(2*dim*n)^2));

J(1:dim,1:dim) = speye(dim);
J(dim+1:2*dim,1+2*(n+1)*dim-dim:2*(n+1)*dim) = speye(dim);

for j = 1:n+1   % block column
    D = [ JFu(t(j),u(:,j),v(:,j)) , JFv(t(j),u(:,j),v(:,j)) ; 
          JGu(t(j),u(:,j),v(:,j)) , JGv(t(j),u(:,j),v(:,j)) ];
    for i = 1:n % block row
        J(2*dim+1+(i-1)*2*dim:2*dim+i*2*dim,1+(j-1)*2*dim:j*2*dim) = C(j,i+1)*D;
    end
end

M12 = blkdiag(M1,M2);
N = -eye(n+1);
N(1,:) = 0;
N(2:end,1) = 1;
J = J + kron(sparse(N),M12);