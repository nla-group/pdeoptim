%% Nonlinear 2D Schloegl control problem with spectral space discr.
% Testing error estimator.
% generates Figure 3.4 in the paper

mydefaults

colorOrder = get(gca, 'ColorOrder'); close all
markers = { ':','-','--*','--s','-d', ':v' };
runs = 0; % nr of runs for timings. can reduce to 0 for just producing plots
nx = 32; % inner grid points in space

%% Initialisation of space discretization
% Set up FD differentiation matrix L on [-1,1]x[-1,1].
T = 2;        % time interval of integration [0,T]
bet = 0.01;   % penalisation parameter
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
JFu = @(t,u,v) L+spdiags(1-3*u.^2,0,nx^2,nx^2);
%JFv = @(t,u,v) jac_v(F,t,u,v);
JFv = @(t,u,v) spdiags(1/bet*ones(nx^2,1),0,nx^2,nx^2);

%JGu = @(t,u,v) jac_u(G,t,u,v);
JGu = @(t,u,v) spdiags(1+6*u.*v,0,nx^2,nx^2);
%JGv = @(t,u,v) jac_v(G,t,u,v);
JGv = @(t,u,v) -L+spdiags(-1+3*u.^2,0,nx^2,nx^2);

%%
n = 5;
tstype = 2;
opts.uex = uex;
opts.vex = vex;
opts.JFu = JFu;
opts.JFv = JFv;
opts.JGu = JGu;
opts.JGv = JGv;
opts.tol = 1e-16;
method = 3;         % iterative solve with FD approx of Jacobian+dyn forcing
iters = 11;         % max number of Newton iterations
opts.verbose = 0;
opts.waitbar = 0;

for j = 1:3  % lag parameter for error estimation
    
    opts.errorest = j; % >=1 will return absolute error estimate
    
    %% CALL TO THE SOLVER
    [t,u,v,out] = solver(F, G, M1, M2, u0, vT, T, n, tstype, method, iters, opts);
    
    figure(102)
    if j == 1
        semilogy(0:length(out.nrmerru)-1,out.nrmerru,markers{4},'Color',colorOrder(4,:));
        hold on
    end
    semilogy(0:length(out.nrmerru_est)-1,out.nrmerru_est,markers{j},'Color','k');
    figure(103)
    if j ==1
        semilogy(0:length(out.nrmerrv)-1,out.nrmerrv,markers{4},'Color',colorOrder(4,:));
        hold on
    end
    semilogy(0:length(out.nrmerrv_est)-1,out.nrmerrv_est,markers{j},'Color','k');
    
end


%%
figure(102)
ylim([1e-6,2])
legend('error','estimate (s=1)','estimate (s=2)','estimate (s=3)','Location','SouthWest')
title('Schloegl FD: error in computed state')
xlabel('Newton iteration')
mypdf('example_Schloegl_FD_32_err_1',.75,1.3)
figure(103)
ylim([1e-6,2])
legend('error','estimate (s=1)','estimate (s=2)','estimate (s=3)','Location','SouthWest')
title('Schloegl FD: error in computed adjoint')
xlabel('Newton iteration')
mypdf('example_Schloegl_FD_32_err_2',.75,1.3)
