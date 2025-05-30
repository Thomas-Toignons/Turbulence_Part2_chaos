close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
dt = 0.1;                      % time step size for time integration
epsilon = 1e-6;                % perturbation amplitude for computing J
tau = 2;                       % time intervals of computing J
N_norm = 5000;                 % number of re-normalizations
N_exp = 10;                    % number of Lyapunov exponents

%% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

%% compute Lyapunov exponents
Q = zeros(length(v0),N_exp);   % allocate matrix of orthonormal vectors
for i = 1:N_exp
    Q(i,i) = 1;                % initialize orthonormal vectors
end

X = zeros(N_exp,N_norm);       % allocate matrix of Lyapunov exp. history
t = zeros(N_norm,1);           % allocate vector of normalization instances


%%
linestyles = {'-',':'};

figure
for i = 1:N_norm
    J = Jacobian(v0,tau,epsilon,dt,L,N,symm); % compute Jacobian matrix
    
    %---To complete---
    V = J*Q; % map deviation vector
    [Q,R] = qr(V, 0); % QR decomposition

    if i==1
        X(:,i) = log(abs(diag(R))); % update X
    else
        X(:,i) = X(:,i-1) + log(abs(diag(R))); % update X
    end
    %---To complete---

    [v0,~] = KSE_integrate(v0,tau,dt,0,L,N,symm); % evolve the orbit and update v0
    t(i) = i*tau;


    
    if(rem(i,20)==0)           % update figure every 20 re-normalizations
        clf; grid on; hold on
        for q = 1:N_exp
            idx = mod(q,2)+1;
            plot(t(1:i),X(q,1:i)./(tau*(1:i)),'LineWidth',2, ...
                Linestyle=linestyles{idx})
        end
        xlabel('Iteration', Interpreter='latex'); ylabel('\chi_i')
        legend(compose('\\chi_{%d}', 1:N_exp))%%
        drawnow
    end
end
exportgraphics(gcf,"../figures/lyapunov.png", Resolution=600)
%%
chi = X(:, N_norm)/(tau*N_norm)