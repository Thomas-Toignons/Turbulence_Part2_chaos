%%% DESCRIPTION -----------------------------------------------------------
%   Jacobian of the KSE flow map: J=dv(t)/dv(0)


%%% INPUTS ----------------------------------------------------------------
%   v0          reference point of the Jacobian (column state vector)
%   t           integration time interval
%   epsilon     perturbation magnitude for finite difference derivatrives
%   dt          step size in time integrations
%   L           domain length
%   N           spatial resolution
%   symm        center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   J           Jacobian matrix


function J = Jacobian(v0,t,epsilon,dt,L,N,symm)
    n = length(v0);
    J = zeros(n,n);
    
    % Loop over each dimension
    for j = 1:n
        e = zeros(n,1);
        e(j) = 1;
        
        % Perturb positively and negatively for central differencing
        v_plus  = v0 + epsilon * e;
        v_minus = v0 - epsilon * e;
        
        % Integrate both perturbed initial conditions
        [v_plus_t,~] = KSE_integrate(v_plus, t, dt, 0, L, N, symm);
        [v_minus_t,~] = KSE_integrate(v_minus, t, dt, 0, L, N, symm);
        
        % Central difference approximation of column j
        J(:, j) = (v_plus_t - v_minus_t) / (2 * epsilon);
    end
end