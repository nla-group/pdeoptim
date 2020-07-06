%% Nonlinear 2D Schloegl control problem with FD space discr.
% this example visualises the eigenvals of the Schur complement S
% and its preconditioned version
% generates Figure 3.1 in the paper

mydefaults

%% Initialisation of space discretization
% Set up FD differentiation matrix L on [-1,1]x[-1,1].
T = 2;         % time interval of integration [0,T]
bet = 0.01;  % penalisation parameter
nx = 16;       % interior grid pts in space
x = linspace(-1,1,nx+2); x = x(2:end-1);
[x1,x2] = meshgrid(x); x1 = x1(:); x2 = x2(:);
L = (nx+1)^2/4*(-2*eye(nx)+diag(ones(nx-1,1),1)+diag(ones(nx-1,1),-1));
I = speye(nx);
L = kron(L,I) + kron(I,L);  % Laplacian on full grid

%% Initial/final values and desired state

a = 2*exp(T)/(bet*(pi^2-2)); b = -2/(bet*pi^2);
constant = 0.1;
u0 = constant*(a+b)*cos(pi*x1/2).*cos(pi*x2/2);
uhat = @(t) constant*(a+(pi^2/2-1)*exp(T)+(2+b-pi^2/2)*exp(t))*cos(pi*x1/2).*cos(pi*x2/2) ...
    + 3*constant^3*(exp(T)-exp(t))*(a+b*exp(t))^2*cos(pi*x1/2).^3.*cos(pi*x2/2).^3;
vT = zeros(nx^2,1);
z = @(t) constant^3*(a+b*exp(t))^3*cos(pi*x1/2).^3.*cos(pi*x2/2).^3;

%% exact solution
uex = @(t) constant*(a+b*exp(t))*cos(pi*x1/2).*cos(pi*x2/2);
vex = @(t) constant*(exp(T)-exp(t)).*cos(pi*x1/2).*cos(pi*x2/2);

%% Define nonlinear ODE system 
%
%   u' = K_1(t,u,v)*u - K_2(t,u,v)*v + f(t,u,v),
%   v' = K_3(t,u,v)*u - K_4(t,u,v)*v + g(t,u,v),
%
% where K_j = K_j(t,u,v).

I = speye(nx^2);
K1 = @(t,u,v) sparse(L) + I - spdiags(u.^2,0,length(u),length(u));
K2 = @(t,u,v) -1/bet*I;
K3 = @(t,u,v) I;
K4 = @(t,u,v) sparse(L) + I - 3*spdiags(u.^2,0,length(u),length(u));
f = @(t,u,v) z(t); % - u.^3;
g = @(t,u,v) -uhat(t); % + 3*u.^2.*v;
M1 = I; M2 = I;

F = @(t,u,v) K1(t,u,v)*u - K2(t,u,v)*v + f(t,u,v);
G = @(t,u,v) K3(t,u,v)*u - K4(t,u,v)*v + g(t,u,v);

% Jacobians
%JFu = @(t,u,v) jac_u(F,t,u,v);
JFu = @(t,u,v) L+diag(1-3*u.^2);
%JFv = @(t,u,v) jac_v(F,t,u,v);
JFv = @(t,u,v) diag(1/bet*ones(nx^2,1));

%JGu = @(t,u,v) jac_u(G,t,u,v);
JGu = @(t,u,v) diag(1+6*u.*v);
%JGv = @(t,u,v) jac_v(G,t,u,v);
JGv = @(t,u,v) -L+diag(-1+3*u.^2);

%%
n = 5;
tstype = 2;
opts.uex = uex;
opts.vex = vex;
opts.JFu = JFu;
opts.JFv = JFv;
opts.JGu = JGu;
opts.JGv = JGv;
opts.tol = 1e-5;
opts.plotevs = 1;

method = 2;        % iterative solve using preconditioned GMRES
iters = 3;         % max number of Newton iterations

%% CALL TO THE SOLVER
%
[t,u,v,out] = solver(F, G, M1, M2, u0, vT, T, n, tstype, method, iters, opts);
%



%%
figure(1), axis([-450,0,-300,300]), shg
h = text(-420,-255,['cond(S)=' num2str(out.condS(1),'%.1e')]); set(h,'FontSize',16)
mypdf('example_Schloegl_FD_eig_1',.75,1.3)
figure(3), axis([-450,0,-300,300]), shg
h = text(-420,-255,['cond(S)=' num2str(out.condS(2),'%.1e')]); set(h,'FontSize',16)
mypdf('example_Schloegl_FD_eig_3',.75,1.3)
figure(5), axis([-450,0,-300,300]), shg
h = text(-420,-255,['cond(S)=' num2str(out.condS(3),'%.1e')]); set(h,'FontSize',16)
mypdf('example_Schloegl_FD_eig_5',.75,1.3)

figure(2), axis([0,1.3,-.3,.301]), shg
h = text(.05,-.24,['cond(P^{-1}S)=' num2str(out.condPS(1),'%.1e')]); set(h,'FontSize',16)
mypdf('example_Schloegl_FD_eig_2',.75,1.3)
figure(4), axis([0,1.3,-.3,.301]), shg
h = text(.05,-.24,['cond(P^{-1}S)=' num2str(out.condPS(2),'%.1e')]); set(h,'FontSize',16)
mypdf('example_Schloegl_FD_eig_4',.75,1.3)
figure(6), axis([0,1.3,-.3,.301]), shg
h = text(.05,-.24,['cond(P^{-1}S)=' num2str(out.condPS(3),'%.1e')]); set(h,'FontSize',16)
mypdf('example_Schloegl_FD_eig_6',.75,1.3)