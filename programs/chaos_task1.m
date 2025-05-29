close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
T_trans = 1000.0;              % transient time period
T_study = 250.0;               % analysis time period
dt = 0.1;                      % time step size for time integration
dt_store = 1.0;                % time intervals of storing a snapshot
epsilon = 1e-2;                % relative ampliture of perturbation

%% Step 1: initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

% transient time integration
[v1000,~] = KSE_integrate(v0,T_trans,dt,0,L,N,symm);

%% Step 2: Advance step v1000 for 250 time units
[v1,t] = KSE_integrate(v1000, T_study, dt, dt_store, L, N, symm);

u1 = zeros([64,250]);
for i=1:size(v1,2)
    u1(:,i) = vector2field(v1(:,i),N,symm);
end

%% Step 3: Perturb v1000 and advance for 250 times units
eta = rand(size(v1000));
r = epsilon * (2*eta - 1) .* v1000;

v_perturbed = v1000 + r;

% Advance v_perturbed for 250 time units
[v2,~] = KSE_integrate(v_perturbed,T_study,dt,dt_store,L,N,symm);

u2 = zeros([64,250]);
for i=1:size(v1,2)
    u2(:,i) = vector2field(v2(:,i),N,symm);
end

%% Step 4: Compute difference
diff = abs(u1 - u2);

%% Plots for steps 2, 3, 4
figure('Name','A');
subplot(1,3,1);
imagesc(u1')
c = colorbar;
title('$u_1(x,t)$', Interpreter='latex');
xlabel('$x$', Interpreter='latex')
ylabel('$t$', Interpreter='latex')
set(gca,'YDir','normal')

subplot(1,3,2);
imagesc(u2')
c = colorbar;
title('$u_2(x,t)$', Interpreter='latex');
xlabel('$x$', Interpreter='latex')
ylabel('$t$', Interpreter='latex')
set(gca,'YDir','normal')

subplot(1,3,3);
imagesc(log10(diff'))
c = colorbar;
title('log10($|u_1 - u_2|)(x,t)$', Interpreter='latex');
xlabel('$x$', Interpreter='latex')
ylabel('$t$', Interpreter='latex')
set(gca,'YDir','normal')

exportgraphics(gcf, '../figures/contour.png',Resolution=600)


%%
% Compute initial pointwise difference
initial_diff = abs(u1(:,1) - u2(:,1));   % shape: [64 x 1]

% Create mask: for each (x,t), check if diff >= e * initial_diff(x)
mask = zeros(size(diff));               % shape: [64 x 250]
for i = 1:size(diff,1)
    mask(i,:) = diff(i,:) >= exp(1) * initial_diff(i);
end

% Plot mask
figure('Name','Lyapunov Time Mask');
imagesc(mask')
set(gca, 'YDir', 'normal')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t$', 'Interpreter', 'latex')
title('$|u_1 - u_2| \geq e \cdot |u_1(x,0) - u_2(x,0)|$', 'Interpreter', 'latex')
colorbar;
