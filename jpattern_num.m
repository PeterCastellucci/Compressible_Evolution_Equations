function S=jpattern_num(H0, M, dxi, zeta, L, N, threshold)
%JPattern_NUM   Computes the sparsity pattern of the Jacobian matrix numerically.
%
%   S = jpattern_num(H0, M, dxi, zeta, L, N)
%
%   This function determines the sparsity pattern of the Jacobian matrix
%   corresponding to the system of evolution equations defined in the function
%   `Equations`. The Jacobian is approximated using finite differences via MATLAB’s
%   `numjac` function, which is useful for implicit ODE solvers that benefit from
%   sparsity exploitation (e.g., ode15s).
%
%   Outputs:
%       S     - Sparse binary matrix representing the Jacobian structure (nonzero pattern)
%
%   The Jacobian is computed numerically by perturbing each variable and observing the
%   change in the time derivative produced by the `Equations` function.
%   Only the structure (pattern of nonzero entries) is retained in the output.
%
%   Note:
%       - This function assumes that the `Equations` function returns a dense Jacobian
%         when perturbed, which is then sparsified using `spones`.
%       - `vectorized='on'` indicates that `Equations` can handle vectorized inputs.
%
% Example usage:
%   S = jpattern_num(H0, M, dxi, zeta, L, N);
%   opts.JPattern = S;  % Used in ODE solver options for efficiency

tbase=0; % Base time for Jacobian evaluation
ybase = H0;


% Compute the corresponding derivative vector
ytbase=Equations(tbase,ybase, M, dxi, zeta, L, N, threshold);
fac=[];% Reused scaling factors (optional output of numjac)
thresh=1e-16;% Threshold for numerical differentiation
vectorized='on';% Enable vectorized evaluation of RHS

% Compute the Jacobian using numerical finite differences
[Jac,fac]=numjac(@(t,H) Equations(t,H, M, dxi, zeta, L, N, threshold),tbase,ybase,ytbase,thresh,fac,vectorized);

% Extract sparsity pattern: convert to sparse binary matrix
S=spones(sparse(Jac));

end
