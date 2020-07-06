function [C,w,maxlen] = quadmat2(a,b,x,d,cache)
% [C,w] = quadmat2(a,b,x,d)
% Returns the collocation matrix C for cumulative integration at
% the points x on [a,b]. The barycentric blending parameter is d.
% If f is a vector of function values f(x), then g = C*f should
% be an approximation to g(x) with g = cumsum(f).
%
% Uses the CUMSUM command of Chebfun. For more details about Chebfun see
%
%   T. A. Driscoll, N. Hale, and L. N. Trefethen, editors, Chebfun Guide,
%   Pafnuty Publications, Oxford, 2014.
%
% For details about QUADMAT2 and rational deferred correction, see
%
%   S. G\"{u}ttel and G. Klein: "Efficient high-order rational integration
%   and deferred correction with equispaced data." Electronic Transactions
%   on Numerical Analysis 41, pp. 443-464, 2014.

persistent CACHE

if nargin < 5
    cache = 0;
end

if cache == 0
    CACHE = [];
end

n = length(x) - 1;

if cache
    for jc = 1:length(CACHE)
        % check for cache hit:
        if (CACHE{jc}.bminusa - (b-a))<10*eps && norm(CACHE{jc}.xminusa - (x-a))<1e-14 && CACHE{jc}.d == d,
            C = CACHE{jc}.C;
            w = CACHE{jc}.w;
            maxlen = CACHE{jc}.maxlen;
            disp('quadmat2: cache hit')
            return
        end
    end
end

w = weights(n,d,x);
C = zeros(n+1,n+1);

maxlen = 0;

if n == 0
    C = 1;
    return
end

for k = 1:ceil((n+1)/2)
    fx = 0*x; fx(k) = 1; % k-th Lagrange function
    c = chebfun(@(xx) bcinterpol(w,x,fx,xx.'),[x(1),x(end)]); % ! overriding adaptive choice for nr of pts
    maxlen = max(maxlen,length(c));
    c = cumsum(c); c = c(:);
    C(:,k) = c(x) - c.vals(1);
end
for k = ceil((n+1)/2)+1:n+1  % uses symmetry of the Lagrange function, see KB, BIT
    C(:,k) = C(end,n+1-k+1) - C(end:-1:1,n+1-k+1);
end

if cache
    jc = length(CACHE) + 1;
    % check for cache hit:
    CACHE{jc}.bminusa = b-a;
    CACHE{jc}.xminusa = x-a;
    CACHE{jc}.d = d;
    CACHE{jc}.C = C;
    CACHE{jc}.w = w;
    CACHE{jc}.maxlen = maxlen;
    disp('quadmat2: cache write')
end
end


%%

function w = weights(n,d,x)
% Compute the barycentric FH weights. This function has been written by
% Georges Klein.

for k=0:n
    ji=max(k-d,0);
    jf=min(k,n-d);
    m=1;
    sumcoeff=[];
    product=1;
    for i=ji:jf
        l=1;
        prodterm=[];
        for j=i:i+d
            if(j==k)
                prodterm(l)=1;
            else
                prodterm(l)=(x(k+1)-x(j+1));
            end
            l=l+1;
        end
        product=1/prod(prodterm);
        sumcoeff(m)=((-1)^(i-1))*product;
        m=m+1;
    end
    [Y,I]=sort(abs(sumcoeff));
    Y(:)=sumcoeff(I(:));
    w(k+1)=sum(Y);
end
w=w(:);
end


%%
function ff=bcinterpol(w,x,f,xx)
% Evaluate barycentric rational interpolant. This function has been written
% by Georges Klein.

[mask,index]=ismember(xx,x);
invmask=(mask==0);
xxx=xx(invmask);

ff=zeros(length(xx),1);
numer=zeros(length(xxx),1)';
denom=zeros(length(xxx),1)';
for i=1:length(x)
    temp=w(i)./(xxx-x(i));
    denom=denom+temp;
    if f(i) == 0
        continue
    end
    numer=numer+temp*f(i);
end

ff(invmask)=numer./denom;
ff(mask)=f(index(mask));
end
