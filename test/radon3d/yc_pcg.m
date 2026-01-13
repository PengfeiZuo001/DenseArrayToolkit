function [m, misfit] = yc_pcg(operator, Param, d, m0, Niter_in, Niter_out, verb)
% YC_PCG - Preconditioned Conjugate Gradient solver for sparsity-promoting inverse problems
%
% This function implements a preconditioned conjugate gradient algorithm
% for solving inverse problems with L1 regularization (sparsity promotion).
% The optimization problem is:
%   min || d - F*m ||_2^2 + μ * || m ||_1
%
% The algorithm uses an iterative reweighting approach to approximate the
% L1 norm with a weighted L2 norm, making it suitable for sparse solutions.
%
% INPUTS:
%   operator - Function handle to forward operator F (e.g., @radon3d_op)
%   Param    - Parameter structure for the forward operator
%   d        - Right-hand side data vector/matrix [nt, nhx, nhy] or [nt, nhx]
%   m0       - Initial model estimation [nt, npx, npy] or [nt, nv]
%   Niter_in - Maximum number of inner iterations (conjugate gradient steps)
%   Niter_out- Number of outer iterations (reweighting steps)
%   verb     - Verbosity flag (1 = display progress, 0 = silent)
%
% OUTPUTS:
%   m        - Estimated model that minimizes the objective function
%   misfit   - History of misfit values across iterations
%
% ALGORITHM FEATURES:
%   - Preconditioned conjugate gradient method
%   - Iterative reweighting for L1 regularization
%   - Sparsity-promoting inversion
%   - Robust convergence for ill-posed inverse problems
%
% REFERENCE:
%   Chen, Y., 2018, Automatic velocity analysis using high-resolution 
%   hyperbolic Radon transform, Geophysics
%
% DEMO:
%   test/test_radon_recon_linear.m
%   test/test_radon_recon_hyper.m
%
% Author: Yangkang Chen (original), MATLAB DenseArrayToolkit (enhancements)
% Date: December 2016 (original), Enhanced for clarity and documentation

u = m0;
P = ones(size(u));
kc = 1;
Misfit = [];
m = u;
for l = 1:Niter_out
    di = feval(operator,P.*u,Param,1);
    r = d-di;
    
    g = feval(operator,r,Param,-1);
    g = g.*P;
    s = g;
    gammam = cgdot(g); %r^T_{k-1}r_{k-1}
    k = 1;
    while  k<=Niter_in
        q = feval(operator,P.*s,Param,1);
        den = cgdot(q);
        alpha = gammam/(den+1.e-8);
        u = u+alpha*s;
        r = r-alpha*q;
        misfit(kc) =  cgdot(r);
        g = feval(operator,r,Param,-1);
        g = g.*P;
        gamma = cgdot(g);%r^T_{k}r_{k}
        beta = gamma/(gammam + 1.e-7);
        gammam = gamma;
        s = g + beta*s;
        if verb
            fprintf('Iteration = %d Misfit=%0.5g\n',k,misfit(k));
        end
        k = k + 1;
        kc = kc + 1;
    end
    
    m = P.*u;
%     P=sqrt(abs(m))+0.00001; %or P=sqrt(abs(m))+0.000001; 
    % m=Pu, P=diag(abs(m0)^0.5);
    
    P = abs(m/max(abs(m(:))))+0.001; %sometimes this one is better, why?
end
return


function [out] = cgdot(in)

% Dot product

temp  =   in.*conj(in);
out = sum(temp(:));

return;
