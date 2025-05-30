close all; clear; clc

set(0, 'DefaultAxesFontSize', 18); % Set default font size for axes
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
    plot(u_guess, 'DisplayName','Initial guess', LineWidth=1.5, LineStyle='--')

    if (params(i).flag) == 1
        plot(u_eq, 'DisplayName','Equilibrium', LineWidth=1.5)
    end
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
clear params

figure('Name','D');
for k=1:6
    params(k).k = k;
    u0 = sin(k* (2*pi*x)/L);
    V_guess = field2vector(u0,N,symm);
    [params(k).V_eq, params(k).flag] = search4EQ(V_guess, T_eqb, dt, L, N, symm);
    
    u_guess = vector2field(V_guess,N,symm);
    u_eq = vector2field(params(k).V_eq,N,symm);

    subplot(3,2,k)
    hold on
    plot(u_guess, 'DisplayName','Initial guess', 'LineStyle','--', linewidth=1.5)

    if (params(k).flag) == 1
        plot(u_eq, 'DisplayName','Equilibrium', LineWidth=1.5)
    end

    title(['$\sin($' num2str(k) '$\cdot \frac{2 \pi x}{L})$'], 'Interpreter','latex')
    xlabel('$x$', 'Interpreter','latex')
    ylabel('$u(x)$', 'Interpreter','latex')
end
legend()
exportgraphics(gcf,'../figures/sines.png', Resolution=600)


[E,P,D] = projection([params.V_eq],N,L,symm);
for i=1:length(snapshots)
    params(i).E = E(i);
    params(i).P = P(i);
    params(i).D = D(i);
end


%% Plot the chaotic attractor
[vv,tt] = KSE_integrate(v0,t_study,dt,dt_store,L,N,symm); % integrated state vector for the chaotic attractor plot, dt_store is lower to have smoother curves
[E,P,D] = projection(vv,N,L,symm);

figure(Name='E');
plot3(E,P,D,'Color',[0,0,0,0.25], DisplayName='Attractor') % Chaotic attractor visualization
hold on

plot3([params.E], [params.P], [params.D],'r', DisplayName='Equilibria from sines',...
    LineStyle='none', Marker='.', MarkerSize=30); % E,D,P visualization

plot3(31.17, 40.32, 40.32, 'g', DisplayName='Equilibria from snapshots', ...
    LineStyle='none', Marker='.', MarkerSize=30); % Add manually the 'missing' equilibrium point



xlabel('Energy $E$', Interpreter='latex');
ylabel('Production $P$', Interpreter='latex');
zlabel('Dissipation $D$', Interpreter='latex');
grid on;
legend();
exportgraphics(gcf, '../figures/attractor.png', 'Resolution',600)

%% computing UPOs
clc; clear
set(0, 'DefaultAxesFontSize', 18); % Set default font size for axes


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

ts = [189, 529, 1113, 2087, 2004];
T_guesses = [37, 58, 74, 85, 109];

for i=1:length(ts)
    t = ts(i);
    T_guess = T_guesses(i);
    
    params(i).t = t;
    params(i).T_guess = T_guess;

    V_guess = V(:,(t+1)); %the index corresponding the the t-value is t+1

    [V_UPO, T_UPO, flags] = search4PO(V_guess,T_guess,dt,L,N,symm);
    
    params(i).flag = flags;
    params(i).T_UPO = T_UPO;
    params(i).V_UPO = V_UPO;
    
end

% Determine stability of computed UPOs

for i=1:length(T_guesses)
    V_p = V(:,i);

    JvT = Jacobian(V_p,params(i).T_UPO,epsilon,dt,L,N,symm);

    params(i).floquet = max(abs(eig(JvT)));

    params(i).tf = params(i).T_UPO / params(i).floquet;
    
end

%% Plot the orbits
[Vc,~] = KSE_integrate(v0,t_study,dt,dt_store,L,N,symm);
[E,P,D] = projection(Vc, N, L, symm);

figure('Name','G')
plot3(E,P,D,'Color',[0,0,0,0.25], DisplayName='Attractor') % Chaotic attractor visualization
hold on;
for i=1:length(ts)
    V = params(i).V_UPO;
    T = params(i).T_UPO;
    
    if params(i).flag == 1
        v = KSE_integrate(V,T,dt,0.1,L,N,symm);
        [E,P,D] = projection(v,N,L,symm);
        plot3(E,P,D, linewidth=2, DisplayName=['UPO ' num2str(i)])
    end

end
xlabel('Energy $E$', Interpreter='latex');
ylabel('Production $P$', Interpreter='latex');
zlabel('Dissipation $D$', Interpreter='latex');
legend()
grid on;
exportgraphics(gcf, '../figures/UPOs.png', Resolution=600);

% close all;

%% Plot each orbit individually
figure('Name','Individual Orbits');
hold on;
k = 1;
for i=1:length(ts)
    V = params(i).V_UPO;
    T = params(i).T_UPO;

    if params(i).flag == 1
        subplot(2,2,k)
        v = KSE_integrate(V,T,dt,0.1,L,N,symm);
        [E,P,D] = projection(v,N,L,symm);
        plot3(E,P,D, linewidth=2)
        
        title(['$T_{UPO} = $' num2str(params(i).T_UPO)], Interpreter="latex")
        xlabel('Energy $E$', Interpreter='latex');
        ylabel('Production $P$', Interpreter='latex');
        zlabel('Dissipation $D$', Interpreter='latex');

        k = k+1;
    end
end
exportgraphics(gcf,'../figures/individual_orbits.png', Resolution=600)