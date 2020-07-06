function [ vec ] = fun_newton_residual_fwd( t, M1, M2, F, G, C, C2, u0, vT, uv)
%FUN_NEWTON_RESIDUAL 
%   Using interlaced forward representation
%   Note: C2 is unused!

dim = size(u0,1); n = length(t)-1;

UV = reshape(uv,2*dim,n+1);
u = UV(1:dim,:); v = UV(dim+1:end,:);

for k = 1:length(t)
    kernelf(:,k) = F(t(k),u(:,k),v(:,k));
    kernelg(:,k) = G(t(k),u(:,k),v(:,k));
end

ru = kernelf*C + M1*(u(:,1)-u); ru(:,1) = u(:,1) - u0;
rv = kernelg*C + M2*(v(:,1)-v); rv(:,1) = v(:,end) - vT;

vec = [ ru; rv ];
vec = vec(:);