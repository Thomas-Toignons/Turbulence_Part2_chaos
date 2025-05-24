close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
t_study = 2500.0;              % analysis time period
dt = 0.1;                      % time step size for time integration
dt_store = 1.0;                % time intervals of storing a snapshot
T_max = 120.0;                 % maximum T for recurrent flow analysis
T_eqb = 10.0;                  % time interval for computing equilibria
epsilon = 1e-6;                % perturbation amplitude for computing J

%% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

%% time integrate the KSE
[V,t_vec] = KSE_integrate(v0,t_study,dt,dt_store,L,N,symm);

%% computing equilibria
snapshots = [2,111,555,666,1705,2222];     % Modify the set of snapshot numbers. Do NOT
                               % use the rand or randi command to select
                               % snapshots from the trajectory; choose six 
                               % arbitrary numbers yourself. 
vb = zeros(size(V,1), length(snapshots));
flg = zeros(length(snapshots));
for i=1:length(snapshots)
    Veq = V(:, snapshots(i));
    [vb(i), flg(i)] = search4EQ(Veq, T_eqb, dt, L, N, symm);
end

%% computing UPOs
%%% to be completed