function [t,u,v,out] = solver(F, G, M1, M2, u0, vT, T, n, tstype, method, iters, opts)
% SOLVER - nonlinear PDE-constrained optimization
%
% This code solves a coupled system of ODEs   
%
%          M1*u' = F(t,u,v), u(0) = u0,
%          M2*v' = G(t,u,v), v(T) = vT.
%
% on the time interval [0,T].
%
% Here,    F(t,u,v) and G(t,u,v) are function handles,
%          M1 and M2 are matrices,
%          u0 and vT are vectors, and
%          T is a scalar > 0.
%
% Further, n is the degree of the time collocation scheme, 
%          tstype is the type of time steps (1 equisdistant, 2 Chebyshev),
%          method is the solution method for the Newton linear system,
%          iters is the number of outer Newton iterations to perform.
% 
% The parameter method can be 
% 
%          method = 1: direct solution of the Newton system using backslash
%          method = 2: preconditioned GMRES using the full matrix Jacobian
%          method = 3: preconditioned GMRES using finite differencing 
%                      (Newton-Krylov, the recommended option)
% 
% The optional parameter structure opts contains 
% function handles to the analytic solution (u(t), v(t)) [if available], 
% and function handles to the Jacobians JFu(t,u,v), JFv, JGu, and JGv.
%
% Outputs: time points t at which solution is collocated,
%          corresponding values of u and v (column-wise), and
%          additional out structure with residual and error information.
%
% The details of this algorithm are described in
%
%          A spectral-in-time Newton--Krylov method 
%          for nonlinear PDE-constrained optimization
%          by Stefan Guettel and John W. Pearson, 2020
%
   
if ~isfield(opts,'errorest')
    opts.errorest = 0;  % 0 = error estimation, >=1 = lag of error est.
end
if ~isfield(opts,'plotevs')
    opts.plotevs = 0;
end
if ~isfield(opts,'verbose')
    opts.verbose = 2;  % 0 = no output, 1 = text output, 2 = text + gmres plot 
end
if ~isfield(opts,'waitbar')
    opts.waitbar = 1;  % 0 = no waitbar, 1 = show waitbar
end

% Setup time integration matrices. There are two choices of time
% collocation points, namely equispaced using barycentric rational
% interpolation (tstype = 1), and polynomial interpolation at Chebyshev
% points in time (tstype = 2).
switch tstype
    case 1
        t = linspace(0,T,n+1); d = min(n,10);
    case 2
        t = chebpts(n+1,[0,T]); d = n;
end
t = t(:).';
[C,w,maxlen] = quadmat2(0,T,t,d,0); C = C.'; % forward time integration matrix
I = eye(n+1); I(n+1,:) = -1; I(n+1,n+1) = 0;
C2 = C*I;                                    % backward time integration matrix
                                             % (only used for checking residual norm)

fun = @(x) fun_newton_residual_fwd(t, M1, M2, F, G, C, NaN, u0, vT, x);

% finite difference approximation of the Jacobian of fun(uv):
Jfun = @(uv) myjac(fun,uv);

% analytic form of the Jacobian of fun(uv):
Jfun2 = @(uv) myjac4(t,M1,M2,opts.JFu,opts.JFv,opts.JGu,opts.JGv,C,NaN,u0,vT,uv);

%% initialize (u,v)
u = repmat(u0,1,n+1);
v = repmat(vT,1,n+1);

if isfield(opts,'uex') && isfield(opts,'vex') % exact solution provided
    for j = 1:length(t)
        U(:,j) = opts.uex(t(j));
        V(:,j) = opts.vex(t(j));
    end
end

%% Perform Newton iters
if opts.waitbar
    wb = waitbar(0,'Newton iters');
end
touter = tic; totalmv = 0; nrmrhs_prev = 1;
for it = 1:iters
    if opts.waitbar
        waitbar(it/iters,wb)
    end
    % compute residuals for u and v
    for k = 1:length(t)
        kernelf(:,k) = F(t(k),u(:,k),v(:,k));
        kernelg(:,k) = G(t(k),u(:,k),v(:,k));
    end
    
    % compute residual for uu
    epsu = kernelf*C;
    epsu = epsu - M1*u + repmat(M1*u0,1,length(t));
    out.nrmresf(it) = max(max(abs(epsu)));
    % compute residual for vv
    epsv = kernelg*C2;
    epsv = repmat(M2*vT,1,length(t)) + epsv - M2*v;
    out.nrmresg(it) = max(max(abs(epsv)));
    
    % compute relative errors if exact solution provided
    if isfield(opts,'uex') && isfield(opts,'vex') 
        out.nrmerru(it) = max(max(abs(u-U)))/max(max(abs(U)));
        out.nrmerrv(it) = max(max(abs(v-V)))/max(max(abs(V)));
    end
    
    % compute relative error estimators using lag parameter opts.errorest
    if opts.errorest
        ind =  rem(it,opts.errorest+1)+1;    % latest index
        ind1 = rem(it+1,opts.errorest+1)+1;  % farthest in past
        HIST_u(:,:,ind) = u;
        HIST_v(:,:,ind) = v;
        if it >= opts.errorest+1
            out.nrmerru_est(it-opts.errorest) = max(max(abs(HIST_u(:,:,ind)-HIST_u(:,:,ind1))))/max(max(abs(HIST_u(:,:,ind))));
            out.nrmerrv_est(it-opts.errorest) = max(max(abs(HIST_v(:,:,ind)-HIST_v(:,:,ind1))))/max(max(abs(HIST_v(:,:,ind))));
        end
    end
    
    % stopping using either error or residual norm 
    if isfield(opts,'tol')
        if isfield(opts,'uex') && isfield(opts,'vex') % exact solution provided
            if max([out.nrmerru(it),out.nrmerrv(it)]) <= opts.tol
                break
            end
        else
            if max([out.nrmresf(it),out.nrmresg(it)]) <= opts.tol
                break
            end
        end
    end
    
    dim = size(u,1);
    uv = [ u ; v ]; uv = uv(:); % interlaced representation
    
    % NEWTON SOLVE
    
    if method == 1  % DIRECT SOLVE
        J = Jfun2(uv);
        rhs = fun(uv);
        step = J\rhs;
    end
    
    if method == 2  % GMRES SOLVE
        
        % Below is the preconditioner of the Schur complement
        % as described in "A preconditioner for the Schur
        % complement using Sherman-Morrison".
        J = Jfun2(uv);
        rhs = fun(uv);
        
        % transform J by adding trailing columns to the first
        K = eye(n+1); K(:,1) = 1; K = kron(K,speye(2*dim));
        JK = J*K;
        
        % Schur complement system matrix (need to solve for this one), 
        % S = D - C*B in paper notation
        S = JK(2*dim+1:end,2*dim+1:end) - JK(2*dim+1:end,1:2*dim)*JK(1:2*dim,2*dim+1:end);
            
        % preconditioner for S
        j = 1;
        Mhat = [ opts.JFu(t(j),u(:,j),v(:,j)) - M1 , opts.JFv(t(j),u(:,j),v(:,j)) ;
                 opts.JGu(t(j),u(:,j),v(:,j)) , opts.JGv(t(j),u(:,j),v(:,j)) - M2 ];
        
        % notation used in the SM formulation
        sC = sum(C); sC = sC(2:end).';
        ctilde = (C(2:end,2:end).')\sC;
        
        %Dhat = kron(C(2:end,2:end)',Mhat);
        %Chat = kron(sC,Mhat);
        %B = JK(1:2*dim,2*dim+1:end);
        %Precond = Dhat - Chat*B;
        
        % the is the inverse of Precond using the SM formula
        %MM = sparse([ zeros(dim,dim), zeros(dim,dim) ; zeros(dim,dim) eye(dim,dim) ]);
        MM = speye(2*dim,2*dim); MM(1:dim,1:dim) = 0;
        MM = kron(ctilde/(1-ctilde(n)),MM);
        FF = speye(2*dim*n);
        FF(:,end+1-2*dim:end) = FF(:,end+1-2*dim:end) + MM;
        %PrecondInv = FF*kron(inv(C(2:end,2:end).'),inv(Mhat));
        
        %FF = 1;
        
        % can implement action of PrecondInv*x efficiently
        PrecondInv_handle = @(x) FF*reshape((Mhat\reshape(x,2*dim,n))/C(2:end,2:end),2*dim*n,1);
        %PrecondInv_handle = @(x) x;
        
        % solution to transformed system (J*K)*sol = rhs is denoted [X;Y]
        % we first find Y using preconditioned GMRES
        rhs2 = rhs(2*dim+1:end) - JK(2*dim+1:end,1:2*dim)*rhs(1:2*dim); % rhs for S-complement
        
        operator = @(x) PrecondInv_handle(S*x);
        tic
        tolres = 1e-3;
        [Y,FLAG,RELRES,ITER,RESVEC] = gmres(operator,PrecondInv_handle(rhs2),[],tolres,100);  %
        totalmv = totalmv + length(RESVEC) - 1;
        
        if opts.verbose >= 1
            disp(['GMRES solve at Newton iter ' num2str(it) ': tolres=' num2str(tolres) ', matvecs=' num2str(length(RESVEC)-1) ', time=' num2str(toc) ])
        end
        if opts.verbose >= 2
            figure(99)
            semilogy(RESVEC), hold on
            title('GMRES convergence on precond Schur complement')
        end
            
        % now find X (easy)
        X = rhs(1:2*dim) - JK(1:2*dim,2*dim+1:end)*Y;
        
        % finally undo column transform
        step = K*[X;Y];  % compare this to J\fun(uv)
        
        if opts.plotevs   % plot evs of original S and preconditioned S
            PrecondInv = FF*kron(inv(C(2:end,2:end).'),inv(Mhat)); % form matrix explicitly
            eeS = eig(full(S));
            eePS = eig(full(PrecondInv*S));
            %eePS = eig(full(S*PrecondInv));
            out.condS(it) = cond(full(S));
            out.condPS(it) = cond(full(PrecondInv*S));
            figure
            plot(eeS,'b.')
            grid on
            title(['evs of Schur compl at Newton iter ' num2str(it)])
            figure
            plot(eePS,'b.')
            grid on
            title(['preconditioned Schur compl at Newton iter ' num2str(it)])
            shg
        end
        
    end
    
    
    if method == 3  % GMRES SOLVE WITHOUT FORMING FULL JACOBIAN
        
        % Below is the preconditioner of the Schur complement
        % as described in "A preconditioner for the Schur
        % complement using Sherman-Morrison".
        
        %J = Jfun2(uv); % no longer needed with finite diff.
        
        rhs = fun(uv);
        
        nrmrhs = norm(rhs);
        % forcing term
        tolres = min(0.1,0.9*nrmrhs^2/nrmrhs_prev^2);
        nrmrhs_prev = nrmrhs;
        
        % transform J by adding trailing columns to the first
        K = eye(n+1); K(:,1) = 1; K = kron(K,speye(2*dim));
        
        %JK = J*K; % no longer needed with finite diff.
        
        % Schur complement system matrix (need to solve for this one)
        %S = JK(2*dim+1:end,2*dim+1:end) - JK(2*dim+1:end,1:2*dim)*JK(1:2*dim,2*dim+1:end);
        
        % get (action of S) by finite differencing of residual function
        
        % recall: S*v = D*v - C*B*v
        
        % just testing to get D*v right
        %vec = [zeros(2*dim,1) ; rhs(2*dim+1:end)]; 
        %h = diffh(uv,vec); dvec = (fun(uv+h*vec) - fun(uv))/h;
        %dvec2 = JK(2*dim+1:end,2*dim+1:end)*rhs(2*dim+1:end);
        %norm(dvec(2*dim+1:end)-dvec2)/norm(dvec2)
        
        % construct handle to D*v
        E = [ sparse(2*n*dim,2*dim) , speye(2*n*dim) ];
        hD = @(x) E*fdiffh(fun,uv,[zeros(2*dim,1) ; x],rhs); % rhs = fun(uv)
        %norm(hD(rhs(2*dim+1:end))-dvec2)/norm(dvec2)
        
        % just testing to get C*v right
        %vec = [randn(2*dim,1) ; zeros(2*n*dim,1) ]; 
        %vec = K*vec;
        %dvec = E*fdiffh(fun,uv,vec);
        %dvec2 = JK(2*dim+1:end,1:2*dim)*vec(1:2*dim);
        %norm(dvec-dvec2)/norm(dvec)
        
        % construct handle to C*v
        hC = @(x) E*fdiffh(fun,uv,K*[x ; zeros(2*n*dim,1)],rhs); % rhs = fun(uv))
        %norm(hC(vec(1:2*dim))-dvec2)/norm(dvec2) 
        
        % construct B (can be done without using JK)
        %B = JK(1:2*dim,2*dim+1:end);
        B = sparse(2*dim,2*n*dim);
        B(dim+1:2*dim,2*dim*n-dim+1:2*dim*n) = speye(dim);
        
        % finally, here is the action of S = D - C*B
        hS = @(x) hD(x) - hC(B*x);
        %vec = randn(2000,1);
        %norm(hS(vec) - S*vec)/norm(S*vec)
        

        % preconditioner for S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = ceil(n/2)+1; % time point in the middle of the interval
        alph = 1.0;
        Mhat = [ opts.JFu(t(j),u(:,j),v(:,j)) - alph*M1 , opts.JFv(t(j),u(:,j),v(:,j)) ;
                 opts.JGu(t(j),u(:,j),v(:,j)) , opts.JGv(t(j),u(:,j),v(:,j)) - alph*M2 ];
         
        % notation used in the SM formulation
        sC = sum(C); sC = sC(2:end).';
        ctilde = (C(2:end,2:end).')\sC;
        
        % the is the inverse of Precond using the SM formula
        %MM = [ sparse(dim,dim), sparse(dim,dim) ; sparse(dim,dim) speye(dim,dim) ];
        MM = speye(2*dim,2*dim); MM(1:dim,1:dim) = 0;
        MM = kron(ctilde/(1-ctilde(n)),MM);
        FF = speye(2*dim*n);
        FF(:,end+1-2*dim:end) = FF(:,end+1-2*dim:end) + MM;
        
        % can implement action of PrecondInv*x efficiently
        dMhat = decomposition(Mhat); 
        PrecondInv_handle = @(x) FF*reshape((dMhat\reshape(x,2*dim,n))/C(2:end,2:end),2*dim*n,1);
        
        %PrecondInv_handle = @(x) x;
        
        % solution to transformed system (J*K)*sol = rhs is denoted [X;Y]
        % we first find Y using preconditioned GMRES
        %rhs2 = rhs(2*dim+1:end) - JK(2*dim+1:end,1:2*dim)*rhs(1:2*dim); % rhs for S-complement
        %rhs2 = rhs(2*dim+1:end) - hC(rhs(1:2*dim)); % get rhs2 by finite diff
        rhs2 = rhs(2*dim+1:end); % if initial/final conditions are already satisfied, rhs(1:2*dim) = 0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        operator = @(x) PrecondInv_handle(hS(x));
        
        tinner = tic;
        [Y,RESVEC] = mygmres(operator,PrecondInv_handle(rhs2),tolres,101);  
        out.gmresit(it) = length(RESVEC)-1;
        totalmv = totalmv + length(RESVEC)-1;
        if opts.verbose >= 1
            disp(['GMRES solve at Newton iter ' num2str(it) ': tolres=' num2str(tolres) ', matvecs=' num2str(length(RESVEC)-1) ', time=' num2str(toc(tinner)) ])
        end
        if opts.verbose >= 2
            figure(99)
            semilogy(RESVEC), hold on
            title('GMRES convergence on precond Schur complement')
        end
        % now find X (easy)
        %X = rhs(1:2*dim) - JK(1:2*dim,2*dim+1:end)*Y;
        %X = rhs(1:2*dim) - B*Y;
        X = -B*Y;   % if initial/final conditions are satisfied, rhs(1:2*dim) = 0
        
        % finally undo column transform
        step = K*[X;Y];  % compare this to J\fun(uv)        
    end

    uv_new = uv - step;
    UV = reshape(uv_new,2*dim,n+1);
    u = UV(1:dim,:); v = UV(dim+1:end,:);
    
end
out.touter = toc(touter);
if opts.verbose >= 1
    disp(['DONE. Total matvec=' num2str(totalmv) ', Total time=' num2str(out.touter) ])
end

if opts.waitbar
    close(wb)
end
end


function d = fdiffh(fun,uv,vec,funuv)
% finite differencing for function fun at uv in the direction vec
% by C. T. Kelley, April 27, 2001

    if nargin < 4
        funuv = fun(uv);
    end
    nrm = norm(vec);
    if nrm==0
        d = zeros(size(vec,1),1);
        return
    end
    h = 1.d-7; xs = (uv'*vec)/nrm;
    if xs ~= 0.d0
         h = h*max(abs(xs),1.d0)*sign(xs);
    end
    h = h/nrm;
    d = (fun(uv+h*vec) - funuv)/h;
end