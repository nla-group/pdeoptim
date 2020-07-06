%% Nonlinear Reaction-diffusion problem
% generates Figure 3.5 and Table 3.3 in the paper

mydefaults
colorOrder = get(gca, 'ColorOrder'); close all
markers = { '-o',':x','-.*','--s','-d', ':v' };
runs = 0; % nr of runs for timings. can reduce to 0 for just producing plots

Nx = [  4, 8, 16 , 32 , 64 , 128]; % inner grid points in space

for j = 1:6
    nx = Nx(j);

    %%
    % Initialisation of space discretization
    % Set up 2d spectral differentiation matrix L on [-1,1]x[-1,1].
    T = 1;       % time interval of integration [0,T]
    x = linspace(-1,1,nx);
    [x1,x2] = meshgrid(x); x1 = x1(:); x2 = x2(:);
    L = (nx-1)^2/4*(-2*eye(nx)+diag(ones(nx-1,1),1)+diag(ones(nx-1,1),-1));
    I = speye(nx);
    L = kron(L,I) + kron(I,L);  % Laplacian on full grid

    % Numeration of nodes:
    %
    % x2 ^
    %    |  4   8  12  16
    %    |  3   7  11  15
    %    |  2   6  10  14
    %    |  1   5   9  13
    %    +-------- x1 ---->

    % restriction to interior grid points  (6 7 10 11)
    RI = eye(nx); RI = RI(2:end-1,:); RI = kron(RI,RI);

    % restriction to boundary points  (1 2 3 4 5 8 9 12 13 14 15 16)
    RB = zeros(2,nx); RB(1,1) = 1; RB(2,nx) = 1;
    RB = blkdiag(eye(nx),kron(eye(nx-2),RB),eye(nx));

    % Normal derivative operator for the grid points
    N = zeros(4*(nx-1),nx^2);

    for i = 1:nx-1
       % Left, working downwards (from top but one)
       N(nx-i,nx-i) = 3/2*((nx-1)/2);
       N(nx-i,2*nx-i) = -2*((nx-1)/2);
       N(nx-i,3*nx-i) = 1/2*((nx-1)/2);

       % Right, working upwards (from bottom but one)
       N(end-nx+i+1,(nx-1)*nx+i+1) = 3/2*((nx-1)/2);
       N(end-nx+i+1,(nx-2)*nx+i+1) = -2*((nx-1)/2);
       N(end-nx+i+1,(nx-3)*nx+i+1) = 1/2*((nx-1)/2);

       % Bottom, working rightwards (from left but one)
       N(nx+2*i-1,i*nx+1) = 3/2*((nx-1)/2);
       N(nx+2*i-1,i*nx+2) = -2*((nx-1)/2);
       N(nx+2*i-1,i*nx+3) = 1/2*((nx-1)/2);

       % Top, working leftwards (from right but one)
       N(3*nx-2-2*i,(nx-i)*nx) = 3/2*((nx-1)/2);
       N(3*nx-2-2*i,(nx-i)*nx-1) = -2*((nx-1)/2);
       N(3*nx-2-2*i,(nx-i)*nx-2) = 1/2*((nx-1)/2);
    end

    % If the following is commented out: at Bottom-left approximate d/dn by
    % -d/dx1, at Top-right approximate by d/dx1; at Bottom-right approximate
    % by -d/dx2; at Top-left approximate by d/dx2

    % Bottom-left (this approximates -1/2*d/dx1-1/2*d/dx2 = constant*d/dn)
    N(1,1) = 3/2*((nx-1)/2);
    N(1,2) = -1*((nx-1)/2);
    N(1,nx+1) = -1*((nx-1)/2);
    N(1,3) = 1/4*((nx-1)/2);
    N(1,2*nx+1) = 1/4*((nx-1)/2);

    % Top-right (this approximates 1/2*d/dx1+1/2*d/dx2 = constant*d/dn)
    N(end,end) = 3/2*((nx-1)/2);
    N(end,end-1) = -1*((nx-1)/2);
    N(end,end-nx) = -1*((nx-1)/2);
    N(end,end-2) = 1/4*((nx-1)/2);
    N(end,end-2*nx) = 1/4*((nx-1)/2);

    % Bottom-right (this approximates 1/2*d/dx1-1/2*d/dx2 = constant*d/dn)
    N(end-nx+1,end-nx+1) = 3/2*((nx-1)/2);
    N(end-nx+1,end-2*nx+1) = -1*((nx-1)/2);
    N(end-nx+1,end-nx+2) = -1*((nx-1)/2);
    N(end-nx+1,end-3*nx+1) = 1/4*((nx-1)/2);
    N(end-nx+1,end-nx+3) = 1/4*((nx-1)/2);

    % Top-left (this approximates -1/2*d/dx1+1/2*d/dx2 = constant*d/dn)
    N(nx,nx) = 3/2*((nx-1)/2);
    N(nx,2*nx) = -1*((nx-1)/2);
    N(nx,nx-1) = -1*((nx-1)/2);
    N(nx,3*nx) = 1/4*((nx-1)/2);
    N(nx,nx-2) = 1/4*((nx-1)/2);

    %% constants and exact solution
    D1 = 0.5; D2 = 1;
    alphay = 1; alphaz = 1; alphau = 0.01;
    k1 = 1; k2 = 1;
    gamma1 = 0.4; gamma2 = 0.6;
    epsilon = 0;

    uex = @(t) [ 1/8*(exp(T)-exp(t))*(x1.^2 + x2.^2);
        1/8*(exp(T)-exp(t))*x1.^0 ];
    vex = @(t) [ alphau*D1/4*(exp(T)-exp(t))*x1.^0 ;
        alphau*D1/4*(exp(T)-exp(t))*(1+cos(pi*x1)).*(1+cos(pi*x2)) ];

    %% problem specification
    yQ = @(t) 1/alphay*(  alphau*D1*exp(t)/4 + (exp(T)-exp(t))*(alphay*(x1.^2+x2.^2)/8+k1*alphau*D1/4) + ...
        alphau*D1/32*(exp(T)-exp(t))^2*(gamma1+gamma2*(1+cos(pi*x1)).*(1+cos(pi*x2)))   );

    zQ = @(t) 1/alphaz*(  alphau*D1/4*exp(t)*(1+cos(pi*x1)).*(1+cos(pi*x2)) + ...
        (exp(T)-exp(t))*(alphaz/8 + pi^2*alphau*D1*D2/4*(2*cos(pi*x1).*cos(pi*x2)+cos(pi*x1)+cos(pi*x2)) + k2*alphau*D1/4*(1+cos(pi*x1)).*(1+cos(pi*x2))) + ...
        alphau*D1/32*(exp(T)-exp(t))^2*(x1.^2+x2.^2).*(gamma1+gamma2*(1+cos(pi*x1)).*(1+cos(pi*x2))));

    y0 = 1/8*(exp(T)-1)*(x1.^2 + x2.^2); % initial condition for y
    z0 = (exp(T)-1)/8*x1.^0; % initial condition for z
    pT = 0*x1;               % final condition for p (reg parameter taken as 0)
    qT = 0*x2;               % final condition for q
    u0 = [ y0 ; z0 ];
    vT = [ pT ; qT ];

    %% Define nonlinear ODE system
    %
    %   M_1 u' = K_1(t,u,v)*u - K_2(t,u,v)*v + f(t,u,v),
    %   M_2 v' = K_3(t,u,v)*u - K_4(t,u,v)*v + g(t,u,v),
    %
    % where K_j = K_j(t,u,v),   u = [y;z],   v = [p;q].

    I = speye(nx^2);

    spzeros = @(m,n) spalloc(m,n,0); 
    N = sparse(N); RI = sparse(RI); RB = sparse(RB);

    K1 = @(t,u,v) [ RI*(D1*L - k1*I - gamma1*spdiags(u(nx^2+1:end),0,nx^2,nx^2)) ,  spzeros(size(RI,1),nx^2)                                   ;
        D1*N                                                         ,  spzeros(size(N,1),nx^2)                                    ;
        spzeros(size(RI,1),nx^2)                                     ,  RI*(D2*L - k2*I - gamma2*spdiags(u(1:nx^2),0,nx^2,nx^2)) ;
        spzeros(size(N,1),nx^2)                                      ,  D2*N+epsilon*RB                                          ];
    K2 = @(t,u,v) [ spzeros(size(RI,1), 2*nx^2) ;
        RB/alphau                 , spzeros(size(RB,1),nx^2) ;
        spzeros(nx^2,2*nx^2)                               ];
    K3 = @(t,u,v) [ RI*(alphay*I)  , spzeros(size(RI,1),nx^2)  ;
        spzeros(size(N,1),2*nx^2)                              ;
        spzeros(size(RI,1),nx^2)   , RI*(alphaz*I)             ;
        spzeros(size(N,1),2*nx^2)                              ];
    K4 = @(t,u,v) [ RI*(D1*L - k1*I - gamma1*spdiags(u(nx^2+1:end),0,nx^2,nx^2)) , -RI*gamma2*spdiags(u(nx^2+1:end),0,nx^2,nx^2)                 ;
        D1*N                                                         ,  spzeros(size(N,1),nx^2)                                    ;
        -RI*gamma1*spdiags(u(1:nx^2),0,nx^2,nx^2)                    ,  RI*(D2*L - k2*I - gamma2*spdiags(u(1:nx^2),0,nx^2,nx^2)) ;
        spzeros(size(N,1),nx^2)                                      ,  D2*N+epsilon*RB                                          ];

    ix1 = RI*x1; ix2 = RI*x2;

    f = @(t,u,v) [
        -1/8*exp(t)*(ix1.^2+ix2.^2)+(exp(T)-exp(t))*(-D1/2+k1/8*(ix1.^2+ix2.^2)) + gamma1/64*(exp(T)-exp(t))^2*(ix1.^2+ix2.^2);
        zeros(size(N,1),1);
        -1/8*exp(t)+(exp(T)-exp(t))*k2/8 + gamma2/64*(exp(T)-exp(t))^2*(ix1.^2+ix2.^2);
        zeros(size(N,1),1) ] ;
    g = @(t,u,v) [
        RI*(-alphay*yQ(t));
        zeros(size(N,1),1);
        RI*(-alphaz*zQ(t));
        zeros(size(N,1),1) ];

    M1 = [ RI , spzeros(size(RI,1),nx^2) ;
        spzeros(size(N,1),2*nx^2)     ;
        spzeros(size(RI,1),nx^2) , RI ;
        spzeros(size(N,1),2*nx^2)     ];
    M2 = M1;

    %%
    F = @(t,u,v) K1(t,u,v)*u - K2(t,u,v)*v + f(t,u,v);
    G = @(t,u,v) K3(t,u,v)*u - K4(t,u,v)*v + g(t,u,v);

    % Jacobians
    JFu = @(t,u,v) jac_u(F,t,u,v);
    JFu = @(t,u,v) [ RI*(D1*L - k1*I - gamma1*spdiags(u(nx^2+1:end),0,nx^2,nx^2)) ,  RI*(- gamma1*spdiags(u(1:nx^2),0,nx^2,nx^2))             ;
        D1*N                                                               ,  spzeros(size(N,1),nx^2)                                    ;
        RI*(- gamma2*spdiags(u(nx^2+1:end),0,nx^2,nx^2))                   ,  RI*(D2*L - k2*I - gamma2*spdiags(u(1:nx^2),0,nx^2,nx^2)) ;
        spzeros(size(N,1),nx^2)                                            ,  D2*N+epsilon*RB                                          ];

    JFv = @(t,u,v) jac_v(F,t,u,v);
    JFv = @(t,u,v) [ spzeros(size(RI,1), 2*nx^2) ;
        -RB/alphau                     , spzeros(size(RB,1),nx^2) ;
        spzeros(nx^2,2*nx^2)                                     ];

    JGu = @(t,u,v) jac_u(G,t,u,v);
    JGu = @(t,u,v) [ RI*(alphay*I)                                            , RI*(gamma1*spdiags(v(1:nx^2),0,nx^2,nx^2) + gamma2*spdiags(v(nx^2+1:end),0,nx^2,nx^2)) ;
        spzeros(size(N,1),2*nx^2)                                                                                                                                         ;
        RI*(gamma1*spdiags(v(1:nx^2),0,nx^2,nx^2) + gamma2*spdiags(v(nx^2+1:end),0,nx^2,nx^2)) , RI*(alphaz*I)                                                          ;
        spzeros(size(N,1),2*nx^2)                                                                                                                                         ];

    JGv = @(t,u,v) jac_v(G,t,u,v);
    JGv = @(t,u,v) [ RI*(-D1*L + k1*I + gamma1*spdiags(u(nx^2+1:end),0,nx^2,nx^2)) ,  RI*gamma2*spdiags(u(nx^2+1:end),0,nx^2,nx^2)              ;
        -D1*N                                                           ,  spzeros(size(N,1),nx^2)                                         ;
        RI*gamma1*spdiags(u(1:nx^2),0,nx^2,nx^2)                        ,  RI*(-D2*L + k2*I + gamma2*spdiags(u(1:nx^2),0,nx^2,nx^2)) ;
        spzeros(size(N,1),nx^2)                                           ,  -D2*N-epsilon*RB                                              ];

    n = 10;
    tstype = 2;
    opts.uex = uex;
    opts.vex = vex;
    opts.JFu = JFu;
    opts.JFv = JFv;
    opts.JGu = JGu;
    opts.JGv = JGv;
    opts.tol = 1e-16;
    opts.plotevs = 0;   % whether to plot eigenvalues of (precond) Jacobian [only with method 2]
    opts.verbose = 0;
    opts.waitbar = 0;
    
    method = 3;         % iterative solve with FD approx of Jacobian+dyn forcing
    iters = 7;          % number of outer Newton iterations

    %% CALL TO THE SOLVER
    [t,u,v,out] = solver(F, G, M1, M2, u0, vT, T, n, tstype, method, iters, opts);
    
    figure(102)
    semilogy(0:length(out.nrmerru)-1,out.nrmerru,markers{j},'Color',colorOrder(j,:));
    hold on
    figure(103)
    semilogy(0:length(out.nrmerrv)-1,out.nrmerrv,markers{j},'Color',colorOrder(j,:));
    hold on
    shg

    outer = 5; % stagnation
    final = max(out.nrmerru(outer), out.nrmerrv(outer));
    total = sum(out.gmresit(1:outer-1));
    inner = sprintf('%d, ', out.gmresit(1:outer-1)); inner = inner(1:end-2);
    
    %% CALL AGAIN TO GET TIMINGS
    iters = outer-1;
    times = [NaN];
    for run = 1:runs
        [t,u,v,out] = solver(F, G, M1, M2, u0, vT, T, n, tstype, method, iters, opts);
        times(run) = out.touter;
    end
    
    disp([ '' , sprintf('%d', nx) , ' & ' , sprintf('%d',outer-1) , ' & ' , inner , ' & ' , sprintf('%d',total) , ' & $' , sprintf('%0.2e', final) , '$ & ' , sprintf('%0.3g', min(times)) , ' \\' ])
    
end

%%
figure(102)
legend('n_x = 4','n_x = 8','n_x = 16','n_x = 32','n_x = 64','n_x = 128','Location','SouthWest')
title('Reactdiff FD: error in computed state')
xlabel('Newton iteration')
mypdf('example_Reactdiff_FD_n_1',.75,1.3)
figure(103)
legend('n_x = 4','n_x = 8','n_x = 16','n_x = 32','n_x = 64','n_x = 128','Location','SouthWest')
title('Reactdiff FD: error in computed adjoint')
xlabel('Newton iteration')
mypdf('example_Reactdiff_FD_n_2',.75,1.3)

%% plotting solutions
ts = 6;
y = u(1:nx^2,ts); z = u(nx^2+1:end,ts); p = v(1:nx^2,ts); q = v(nx^2+1:end,ts);
Y = reshape(y,nx,nx); Z = reshape(z,nx,nx); P = reshape(p,nx,nx); Q = reshape(q,nx,nx);
X1 = reshape(x1,nx,nx); X2 = reshape(x2,nx,nx);
figure(104)
subplot(1,2,1), h = surf(X1,X2,Y,'LineStyle','none'); title('Y'), %zlim([-1,1])
shading interp; lightangle(-85,40); h.FaceLighting = 'gouraud'; h.AmbientStrength = 0.3; h.DiffuseStrength = 0.8; h.SpecularStrength = 0.9; h.SpecularExponent = 25; h.BackFaceLighting = 'unlit';

%subplot(2,2,2), h = surf(X1,X2,Z,'LineStyle','none'); title('Z'), %zlim([-1,1])
%shading interp; lightangle(-85,40); h.FaceLighting = 'gouraud'; h.AmbientStrength = 0.3; h.DiffuseStrength = 0.8; h.SpecularStrength = 0.9; h.SpecularExponent = 25; h.BackFaceLighting = 'unlit';

%subplot(2,2,3), h = surf(X1,X2,P,'LineStyle','none'); title('P'), %zlim([-1,1])
%shading interp; lightangle(-85,40); h.FaceLighting = 'gouraud'; h.AmbientStrength = 0.3; h.DiffuseStrength = 0.8; h.SpecularStrength = 0.9; h.SpecularExponent = 25; h.BackFaceLighting = 'unlit';

subplot(1,2,2), h = surf(X1,X2,Q,'LineStyle','none'); title('Q'), %zlim([-1,1])
shading interp; lightangle(-85,40); h.FaceLighting = 'gouraud'; h.AmbientStrength = 0.3; h.DiffuseStrength = 0.8; h.SpecularStrength = 0.9; h.SpecularExponent = 25; h.BackFaceLighting = 'unlit';


mypdf('example_Reactdiff_FD_n_3',.5,2.2)

