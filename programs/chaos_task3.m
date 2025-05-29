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

%% 4.1 Equilibrium solutions
% STEP 1: RANDOM INITIAL GUESSES

snapshots = [111,256,555,666,1705,2222];

for i=1:6
    params(i).t = snapshots(i)-1;
end

figure('Name','C');
for i=1:length(snapshots)
    t_idx = snapshots(i);
    V_guess = V(:, t_idx);

    [params(i).V_eq, params(i).flag] = search4EQ(V_guess, T_eqb, dt, L, N, symm);
    
    u_guess = vector2field(V_guess,N,symm);
    u_eq = vector2field(params(i).V_eq,N,symm);
    
    subplot(3,2,i)
    hold on
    plot(u_guess, 'DisplayName','Initial guess', linewidth=1.5)
    plot(u_eq, 'DisplayName','Equilibrium', 'LineStyle','--', LineWidth=1.5)
    title(['$t=$' num2str(params(i).t) ' s'], 'Interpreter','latex')
    xlabel('$x$', 'Interpreter','latex')
    ylabel('$u(x)$', 'Interpreter','latex')
end
legend()
exportgraphics(gcf, '../figures/snapshots.png', 'Resolution',600)

[E,P,D] = projection([params.V_eq],N,L,symm);

for i=1:length(snapshots)
    params(i).E = E(i);
    params(i).P = P(i);
    params(i).D = D(i);
end


%% SIN INITIAL GUESSES
flg = zeros(1, 6);
V_eq = zeros(21, 6);

figure('Name','D');
for k=1:6
    u0 = sin(k* (2*pi*x)/L);
    V_guess = field2vector(u0,N,symm);
    [V_eq(:,k), flg(k)] = search4EQ(V_guess, T_eqb, dt, L, N, symm);


    u_guess = vector2field(V_guess,N,symm);
    u_eq = vector2field(V_eq(:,k),N,symm);

    subplot(3,2,k)
    hold on
    plot(u_guess, 'DisplayName','Initial guess', linewidth=1.5)

    if flg(k)~=0
        plot(u_eq, 'DisplayName','Equilibrium', 'LineStyle','--', LineWidth=1.5)
    end

    title(['$\sin($' num2str(k) '$\cdot \frac{2 \pi x}{L})$'], 'Interpreter','latex')
    xlabel('$x$', 'Interpreter','latex')
    ylabel('$u(x)$', 'Interpreter','latex')
end
legend()
exportgraphics(gcf,'../figures/sines.png', Resolution=600)

%% Plot the chaotic attractor
[vv,tt] = KSE_integrate(v0,t_study,dt,0.1,L,N,symm); % integrated state vector for the chaotic attractor plot, dt_store is lower to have smoother curves
[E,P,D] = projection(vv,N,L,symm);

figure(Name='E');
plot3(E,P,D,'Color',[0,0,0,0.25]) % Chaotic attractor visualization

hold on;

[E,P,D]=projection(V_eq,N,L,symm);
plot3(E, P, D,'r', LineStyle='none', Marker='.', MarkerSize=20); % E,D,P visualization

xlabel('Energy $E$', Interpreter='latex');
ylabel('Production $P$', Interpreter='latex');
zlabel('Dissipation $D$', Interpreter='latex');
title('Chaotic attractor of the KSE');
grid on;
exportgraphics(gcf, '../figures/attractor.png', 'Resolution',600)

%% computing UPOs
clc; clear

% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
t_study = 2500.0;              % analysis time period
dt = 0.1;                      % time step size for time integration
dt_store = 1;                % time intervals of storing a snapshot
T_max = 120.0;                 % maximum T for recurrent flow analysis
T_eqb = 10.0;                  % time interval for computing equilibria
epsilon = 1e-6;                % perturbation amplitude for computing J

% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector


% Recurrence indicator
T_tot = 120;
t_tot = t_study;
r = zeros(t_tot, T_tot);


[V,~] = KSE_integrate(v0,t_tot+T_tot,dt,dt_store,L,N,symm);


for t=1:t_tot
    for T=1:T_tot
    idx = t+T;

    V_t = V(:,t);
    V_T = V(:,idx);

    r(t,T) = norm(V_T-V_t)/norm(V_t);

    end
end


figure('Name','F')
imagesc(1:t_tot,1:T_tot,log10(r'))
c = colorbar;
c.Label.String = 'log10($r(t,T)$)';
c.Label.Interpreter = 'latex';
xlabel('$t$', 'Interpreter','latex')
ylabel('$T$', Interpreter='latex')
exportgraphics(gcf,'../figures/recurrence.png', Resolution=600)


%% Search for UPOs
for i=1:5
    params(i).name = i;
end

%%


ts = [189, 529, 1113, 2087, 1366];
T_guesses = [37, 58, 74, 85, 105];
V_UPO = zeros(21, length(ts));
T_UPO = zeros(1, length(ts));
flags = zeros(1, length(ts));


for i=1:length(ts)
    t = ts(i);
    T_guess = T_guesses(i);
    
    params(i).t = t;
    params(i).T_guess = T_guess;

    V_guess = V(:,(t+1)); %the index corresponding the the t-value is t+1

    [V_UPO(:,i), T_UPO(i), flags(i)] = search4PO(V_guess,T_guess,dt,L,N,symm);
    
    params(i).flag = flags(i);
    params(i).T_UPO = T_UPO(i);
    params(i).V_UPO = V_UPO(:,i);
    
end

% Determine stability of computed UPOs
floquets = zeros(1, length(ts));

for i=1:length(T_guesses)
    V_p = V(:,i);
    T_conv = T_UPO(i);
    JvT = Jacobian(V_p,T_conv,epsilon,dt,L,N,symm);
    floquets(i) = max(abs(eig(JvT)));

    params(i).floquet = floquets(i);

    params(i).tf = params(i).T_UPO / params(i).floquet;
    
end

%% Plot the orbits
[Vc,~] = KSE_integrate(v0,t_study,dt,dt_store,L,N,symm);
[E,D,P] = projection(Vc, N, L, symm);

figure('Name','G')
plot3(E,D,P,'Color',[0,0,0,0.25]) % Chaotic attractor visualization
hold on;
for i=1:length(ts)
    V = params(i).V_UPO;
    T = params(i).T_UPO;

    v = KSE_integrate(V,T,dt,0.1,L,N,symm);

    [E,D,P] = projection(v,N,L,symm);
    plot3(E,D,P, linewidth=2)

end
xlabel('Energy $E$', Interpreter='latex');
ylabel('Production $P$', Interpreter='latex');
zlabel('Dissipation $D$', Interpreter='latex');
legend({'KSE attractor', 'UPO 1', 'UPO 2', 'UPO 3', 'UPO 4', 'UPO 5'})
grid on;
exportgraphics(gcf, '../figures/UPOs.png', Resolution=600);

close all;