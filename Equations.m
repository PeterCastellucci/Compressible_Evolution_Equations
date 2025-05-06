% This function implements the discretized equations and boundary conditions
% for time-stepping the system. The domain is split into regions based on Xl:
%
% - When Xl > 0: There are two regions.
%     Region I:   0 < x < Xl — solve for pressure P.
%     Region II:  Xl < x < Xu — solve for B = (1 - F) * P and for F.
%
% - When Xl = 0: The domain is a single region, 0 < x < Xu.
%     In this case, only B and F are evolved.
% 
% This solver supports a time-dependent (nonsteady) source term.
% To use it, define Q as a general function of time t instead of a constant.

function fval = Equations(t,H, M, dxi, zeta, L, N, threshold)

global F0

%Set the source function Q(Omega t)
Q = 1;
% Initialize time derivatives (dHdt) vector
% Last two entries correspond to Xl (lower) and Xu (upper) contact line positions
dHdt = zeros(3*N+2,1);

% Extract contact line positions from state vector
Xu = H(end);
Xl = H(end-1);
dX = Xu - Xl;

% Handle case when lower contact line reaches origin
if Xl < threshold

    % Recover pressure field in Region I: P = B / (1 - F)
    H(2:N-1) = transpose(H(N+2:2*N-1)./(1-H(2*N+2:3*N-1)));

    %Calculate pressure and interface height at origin using boundary
    %condition zeta*(1-F)*P*Px = -Q
    A = 4/3 * H(2) -1/3*H(3) + 4/3*H(2*N+2) -1/3*H(2*N+3);
    a = -3;
    b = 4*H(2) - H(3) -3*(1-A);
    c = (4*H(2) - H(3))*(1-A);
    d = 2*Xu*dxi*Q/zeta;
    
    R = roots([a b c d]);

    diffs = [abs(R - H(2))];
    [~,idx] = min(diffs);

    H(1) = R(idx);% Set origin pressure
    
    %Compute F at the origin
    H(2*N+1) = -H(1) + A;

    %Update global value of F at origin
    F0 = H(2*N+1);
end

%Case where Xl is away from origin (general case)
if Xl >= threshold

    H(2*N+1) = 0; % Enforce F = 0 at lower contact line
    H(3*N) = 1; % Enforce F = 1 at upper contact line
    
    % Recover pressure from B: P = B / (1 - F)
    P = H(N+2:2*N-1)./(1-H(2*N+2:3*N-1));

    %Calculate P at origin from the boundary condition -zeta*P*Px = Q
    a = 3;
    b = -(4*H(2) - H(3));
    c = -2*Xl*dxi*Q/zeta;
    
    H(1) = (-b + sqrt(b^2 -4*a*c))/(2*a);
   
    %Enforce continuity of pressure at lower contact line
    H(N) = 1/(3*(dX/Xl + 1)) *(4*P(1) - P(2) - dX/Xl * (H(N-2) -4*H(N-1)));
    H(N+1) = H(N);
    
    %Calculate pressure at upper contact line (Xu) using boundary condition
    % P + (Px + Fx)*(L-Xu) = 1/zeta - 1
    Pu = 1/(1+3*(L-Xu)/(2*dX*dxi)) *(1/zeta - 1 - (L-Xu)/(2*dX*dxi)*(3 - 4*H(3*N-1) +H(3*N-2) -4*P(end) +P(end-1)));
    P = [H(N), P(:)', Pu];% Extend P vector to include boundary values

    % Compute dXl/dt
    Xlt = - M/(2*dX*dxi) * (-3*P(1) + 4*P(2) - P(3) + 4*H(2*N+2) - H(2*N+3) );
    % Compute dXu/dt
    Xut = - (3*P(end) - 4*P(end-1) + P(end-2))/(2*dxi*dX);
        
    % Store contact line velocities
    dHdt(end) = Xut;
    dHdt(end - 1) = Xlt;
    
    % Time derivatives for pressure field in Region I using central differencing
    for i = 2:N-1
        dHdt(i) = (i-1)*dxi*Xlt*(H(i+1) - H(i-1))/(2*Xl*dxi) + 1/(2*Xl^2*dxi^2) * ((H(i) + H(i+1)) * (H(i+1) - H(i)) - (H(i) + H(i-1)) * (H(i) - H(i-1)));
    end
    
    % Time derivatives for B = (1-F)*P field in Region II using central differencing
    for i = N+2:2*N-1
        q1 = (H(i) + H(i+1))/2 * (P(i+1-N) - P(i-N))/dxi;
    
        q2 = (H(i) + H(i-1))/2 * (P(i-N) - P(i-1-N))/dxi;

        dHdt(i) = ((i-N-1) * dxi * (Xut - Xlt) +Xlt)  * (H(i+1) - H(i-1))/(2*dX*dxi) + 1/(dxi*dX^2) * (q1 - q2);
    end
    
    % Time derivatives for F field in Region II using central differencing
    for i = 2*N+2:3*N-1
        dHdt(i) = ((i-2*N-1) * dxi * (Xut -Xlt) +Xlt) * (H(i+1) - H(i-1))/(2*dxi*dX) + M /(2*dxi^2*dX^2) * ( (H(i)+H(i+1))*(H(i+1) - H(i) + P(i+1-2*N) - P(i-2*N)) - (H(i)+H(i-1))*(H(i) - H(i-1) + P(i-2*N) - P(i-1-2*N)) );
    end

% Special case when lower contact line reaches origin (Xl = 0)  
else
    % Fix Xl = 0 explicitly
    H(end-1) = 0;

    % Recover pressure field: P = B / (1 - F)
    H(2:N-1) = transpose(H(N+2:2*N-1)./(1-H(2*N+2:3*N-1)));

    % Compute value of B at origin using B = (1 - F0)*P0
    H(N+1) = (1-H(2*N+1))*H(1);

    %Calculate pressure at upper contact line (Xu) using boundary condition
    % P + (Px + Fx)*(L-Xu) = 1/zeta - 1
    H(N) = 1/(1+3*(L-Xu)/(2*Xu*dxi)) *(1/zeta - 1 - (L-Xu)/(2*Xu*dxi)*(3 - 4*H(3*N-1) +H(3*N-2) -4*H(N-1) +H(N-2)));

    % Compute dXu/dt
    Xut = - (3*H(N) - 4*H(N-1) + H(N-2))/(2*dxi*Xu);
    dHdt(end) = Xut;
    
    % Time derivatives for B = (1-F)*P field using central differencing
    for i = N+2:2*N-1
        q1 = (H(i) + H(i+1))/2 * (H(i+1-N) - H(i-N))/dxi;
    
        q2 = (H(i) + H(i-1))/2 * (H(i-N) - H(i-1-N))/dxi;

        dHdt(i) = ((i-N-1) * dxi * Xut)  * (H(i+1) - H(i-1))/(2*Xu*dxi) + 1/(dxi*Xu^2) * (q1 - q2);
    end

    % Time derivatives for F field in Region II using central differencing
    for i = 2*N+2:3*N-1
        dHdt(i) = (i-2*N-1) * dxi *Xut * (H(i+1) - H(i-1))/(2*dxi*Xu) + M /(2*dxi^2*Xu^2) * ( (H(i)+H(i+1))*(H(i+1) - H(i) + H(i+1-2*N) - H(i-2*N)) - (H(i)+H(i-1))*(H(i) - H(i-1) + H(i-2*N) - H(i-1-2*N)) );
    end
end
fval = dHdt;
end