% This script numerically solves the compressible evolution equations using the
% method-of-lines (MOL) approach with second-order central differencing in space.
%
% The system evolves two primary quantities:
%   - F: the interface height
%   - B: the depth-integrated density, defined as B = (1 - F) * P, where P is pressure
%
% The spatial domain is divided into two dynamic regions when the lower contact line Xl > 0:
%   - Region I:  0 < x < Xl  → pressure P is evolved here
%   - Region II: Xl < x < Xu → B and F are evolved here
% When Xl = 0, only Region II exists (0 < x < Xu), and both B and F are evolved in that region.
%
% The governing equations and boundary conditions are encoded in the function `Equations`,
% which returns the time derivative of the full system state H, including:
%   - Interface height (F), pressure (P), and depth-integrated density (B)
%   - Dynamic boundary conditions at the moving contact lines (Xl and Xu)
%
% The system can accommodate a time-dependent source term Q(t). To implement this,
% define Q as a function of time within the `Equations` function.
%
% Time integration is performed using a suitable ODE solver (e.g., ode15s),
% applied to the right-hand side provided by the `Equations` function.
%
% Inputs into solver:
%       H0    - Initial state vector (includes F, B, P, Xl, Xu)
%       M     - viscosity ratio parameter
%       dxi   - Spatial grid spacing
%       zeta  - Compressibility parameter
%       L     - Domain length
%       N     - Number of spatial grid points per region

global F0
%Parameters
M = 0.01;
zeta = 1;
L = 100;

% Set the threshold for when Xl reaches the origin (i.e Xl = 0 when Xl <
% threhold)
threshold = 0.001;

% Time discretization step
dt = 0.01;
% Spatial grid spacing in the xi-coordinate (moving frame)
dxi = 0.001;
%Final time
tf = 200;

% Output filename to store simulation data
File1 = 'Test.txt';

% Time domain for integration
tspan = 0:dt:tf;

% Uniform grid in the reference coordinate xi
xi = 0:dxi:1;
% Number of grid points
N = length(xi);

% Initialize state vector: [P, B, F, Xl, Xu]
H0 = zeros(1,3*N+2);

% Initial locations of the lower and upper contact lines
Xl0 = 1;
Xu0 = 3;
dX0 = Xu0 - Xl0;

% Initial condition for F: linear profile between Xl0 and Xu0
H0(2*N+1:3*N) = (Xu0 - Xl0) * xi/dX0;

% Initial condition for pressure P: linear profile which is consistent with
% boundary conditions.
a = 1;
b = (L-Xu0)/dX0 + 1 -1/zeta;
c = -L/zeta;

Pc = (-b + sqrt(b^2 - 4*a*c))/(2*a);

H0(1:N) = Pc - 1/(Pc * zeta) * Xl0 * xi;
P2 = Pc - 1/(Pc * zeta) * (Xl0 + (Xu0-Xl0)* xi);

% Compute initial depth-integrated density B = (1 - F) * P
H0(N+1:2*N) = (1-H0(2*N+1:3*N)).*P2;

% Set initial contact line positions in the state vector
H0(end - 1) = Xl0;
H0(end) = Xu0;

% Compute sparsity pattern for the Jacobian for efficiency in ODE solver
S = jpattern_num(H0, M, dxi, zeta, L, N, threshold);

% Track value of F at the origin to trigger event switching when interface touches x = 0
F0 = 0;

% Set the options for the ODE solver
options = odeset('OutputFcn',@odeprog, 'reltol', 1e-8, 'abstol', 1e-8, "JPattern",S, 'events', @(t,y) events(t,y,L,threshold));

%Integrate the discretized equations, integration stops when Xu reaches the
%outlet or Xl has crossed the origin twice (i.e it has returned to the lower boundary).
[tSol,HSol1, te, ye, ie] = ode15s(@(t,H) Equations(t,H, M, dxi, zeta, L, N, threshold), tspan, H0, options);
rows = size(HSol1,1);

% Display event triggers and time of occurrence
for i = 1:length(ie)
    switch ie(i)
        case 1
            fprintf('Event 1 triggered: Xu exceeded L at t = %.4f\n', te(i));
        case 2
            fprintf('Event 2 triggered: Xl crossed zero at t = %.4f\n', te(i));
        case 3
            fprintf('Event 2 triggered: F0 crossed zero at t = %.4f\n', te(i));
        otherwise
            fprintf('Unknown event triggered at t = %.4f\n', te(i));
    end
end

% Restart integration if interface has crossed the origin twice

if ismember(3,ie)
    
    %Q must match value in equations, e.g Q(Omega t) = t, then put Q =
    %tSol(end)
    Q = 1;
    H0 = HSol1(end,:);
    Xu = H0(end);
    % Compute pressure field from B using the definition B = (1 - F) * P
    H0(2:N-1) = H0(N+2:2*N-1)./(1-H0(2*N+2:3*N-1));
    
    % Use cubic equation root-finding to determine pressure at origin
    A = 4/3 * H0(2) -1/3*H0(3) + 4/3*H0(2*N+2) -1/3*H0(2*N+3);
    a = -3;
    b = 4*H0(2) - H0(3) -3*(1-A);
    c = (4*H0(2) - H0(3)).*(1-A);
    d = 2*Xu*dxi*Q/zeta;
    
    R = roots([a b c d]);

    diffs = abs(H0(2) - R);
    [~,idx] = min(diffs);

    %Origin pressure
    H0(1) = R(idx);
    
    % Compute pressure at the upper contact line from boundary condition
    H0(N) = 1/(1+3*(L-Xu)/(2*Xu*dxi)) *(1/zeta - 1 - (L-Xu)/(2*Xu*dxi).*(3 - 4*H0(3*N-1) +H0(3*N-2) -4*H0(N-1) +H0(N-2)));
    
    % Reset interface height at origin to zero
    H0(2*N+1) = 0;
    P0 = H0(1:N);

    x = Xu*xi;
    
    Xl = 0.01;
    ind = find(x > Xl, 1, 'first');
    H0(end-1) = Xl;

    x1 = Xl*xi;
    x2 = Xl + (Xu-Xl)*xi;
    
    % Interpolate pressure profiles for the one region back to two regions
    P1Interp = interp1(x(1:ind), H0(1:ind), x1, 'spline');
    P2Interp = interp1(x(ind:end), H0(ind:N), x2, 'spline');
    
    % Reconstruct initial condition and time span for second phase of simulation
    H0(1:N) = P1Interp;
    H0(N+1:2*N) = (1-H0(2*N+1:3*N)).*P2Interp;
    tspan2 = rows*dt:dt:tspan(end);
    
    % Solve PDE system again with updated initial conditions
    [tSol2,HSol2, te2, ye2, ie2] = ode15s(@(t,H) Equations(t,H, M, dxi, zeta, L, N, threshold), tspan2, H0, options);
    
    % Report events triggered in the second phase
    for i = 1:length(ie2)
        switch ie2(i)
            case 1
                fprintf('Event 1 triggered: Xu exceeded L at t = %.4f\n', te2(i));
            case 2
                fprintf('Event 2 triggered: Xl crossed zero at t = %.4f\n', te2(i));
            case 3
                fprintf('Event 2 triggered: F0 crossed zero at t = %.4f\n', te2(i));
            otherwise
                fprintf('Unknown event triggered at t = %.4f\n', te2(i));
        end
    end
end

% Combine results from both simulation phases
if ismember(3,ie)
    HSol = vertcat(HSol1, HSol2(2:end,:));
    tSol = [tSol1, tSol2];
else
    HSol = HSol1;
end

% Loop over solution snapshots to apply boundary conditions and reconstruct pressure
for i = 1:size(HSol,1)
    Xu = HSol(i,end);
    Xl = HSol(i,end-1);
    dX = Xu - Xl;
    
    %Q must match value in equations, e.g Q(Omega t) = t, then put Q =
    %tSol(i)
    Q = 1;

    % Case: Xl > 0 — origin not yet reached, standard boundary conditions
    if Xl >= 0.001
        
        HSol(i,2*N+1) = 0;
        HSol(i,3*N) = 1;
        % Swap B field with pressure using B = (1 - F) * P
        HSol(i,N+2:2*N-1) = HSol(i,N+2:2*N-1)./(1-HSol(i,2*N+2:3*N-1));
        
        %Calculate pressure at origin
        a = 3;
        b = -(4*HSol(i,2) - HSol(i,3));
        c = -2*Xl*dxi*Q/zeta;
        
        HSol(i,1) = (-b + sqrt(b.^2 -4*a*c))./(2*a);
        
        % Enforce continuity of pressure at the lower contact line
        HSol(i,N) = 1./(3*(dX./Xl + 1)) .*(4*HSol(i,N+2) - HSol(i,N+3) - dX./Xl .* (HSol(i,N-2) -4*HSol(i,N-1)));
        HSol(i,N+1) = HSol(i,N);
        
        %Upper contact line pressure
        HSol(i,2*N) = 1./(1+3*(L-Xu)./(2*dX*dxi)) .*(1/zeta - 1 - (L-Xu)./(2*dX*dxi).*(3 - 4*HSol(i,3*N-1) +HSol(i,3*N-2) -4*HSol(i,2*N-1) +HSol(i,2*N-2)));

    % Case: Xl = 0 — interface has touched the origin; apply new boundary
    % condition
    else        
        HSol(i,end-1) = 0;
        
        % Reconstruct pressure from B at internal points
        HSol(i,2:N-1) = HSol(i,N+2:2*N-1)./(1-HSol(i,2*N+2:3*N-1));
        
        % Compute interface height and B at origin
        A = 4/3 * HSol(i,2) -1/3*HSol(i,3) + 4/3*HSol(i,2*N+2) -1/3*HSol(i,2*N+3);
        a = -3;
        b = 4*HSol(i,2) - HSol(i,3) -3*(1-A);
        c = (4*HSol(i,2) - HSol(i,3)).*(1-A);
        d = 2*Xu*dxi*Q/zeta;
        
        R = roots([a b c d]);

        diffs = abs(HSol(i,2) - R);
        [~,idx] = min(diffs);
        %Origin pressure
        HSol(i,1) = R(idx);
        
        %Origin F
        HSol(i,2*N+1) = -HSol(i,1) + A;
        
        %Origin B
        HSol(i,N+1) = (1-HSol(i,2*N+1))*HSol(i,1);
        
        %Upper contact line pressure
        HSol(i,N) = 1/(1+3*(L-Xu)/(2*Xu*dxi)) *(1/zeta - 1 - (L-Xu)/(2*Xu*dxi).*(3 - 4*HSol(i,3*N-1) +HSol(i,3*N-2) -4*HSol(i,N-1) +HSol(i,N-2)));
    end
end

% Write final solution matrix to file
writematrix(HSol, File1,'Delimiter',' ')

function [value, isterminal, direction] = events(t,y, L, threshold)
    
    % Event function: returns values and directions for stopping criteria
    % Event 1: Xu reaches outlet (termination)
    % Event 2: Xl hits origin (used for logic flow but doesn't stop integration)
    % Event 3: F at origin becomes zero (termination)
    global F0;
    %Check if Xu has reached outlet
    event1 = y(end) - L;

    %check if Xl has reached origin
    event2 = y(end-1)-threshold;
    
    %check if F at origin has crossed zero
    event3 = F0 - threshold;
    value = [event1; event2; event3];
    isterminal = [1;0;1];

    direction = [1;0;-1];
end
