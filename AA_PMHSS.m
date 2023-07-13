function [xout, k, inner_iter, resvec, relresidual] = AA_PMHSS(C,x0,b,rel_tol,k_max,M0,M1,inner,max_inner_iter,timings)
    %%AA-PHMSS Anderson Accelerated Preconditioned Modified Skew-Herimitian
    %%Splitting 
    %   x = AA_PMHSS(C,x0,b,rel_tol,k_max) attempts to solve the complex
    %   system of linear equations C*x = b
    %   for x.  The N-by-N coefficient matrix C must be complex symmetric
    %   or sym(imag(C)) must be positive semi-definite and the right
    %   hand side column vector b must have length N. 
    %  
    %   AA    (Anderson Acceleration)   
    %           maxit = k_max (argument)
    %           m_k = k_max   (argument)
    %           beta = 0.5 
    %   PMHSS (Preconditioned Modified Skew-Herimitian Splitting) 
    %           M: alpha = 1
    %           P: V = A          
    %   inner solver 
    %           Type = CG 
    %           Tol = 1e-15
    %           maxit = 400
    
    A = real(C);
    B = imag(C);
    V = A;
    alpha = complex(1);
    lhs = 2*(A+B);
    
    
    normb = norm(b);
    m = k_max;
    x = [x0, PMHSS_step(lhs,A,B,b,alpha,V,x0,M0,M1,inner,[],max_inner_iter, x0)]; % Vector of iterates x.
    g = x(:,2) - x0; % Vector of residuals.
    
    X_k = x(:,2) - x(:,1); % Matrix of increments in x.
    G_k = g(:,1); % Matrix of increments in residuals.
    
    residual = norm(G_k(:,1));
    k = 1;
    resvec = [];
    inner_iter = [];
    fx = x(:,2);
    while k <= k_max && residual > rel_tol*normb
        m_k = min(k, m);
     
        gamma_k = G_k\g(:,k);
        x(:,k + 1) = x(:,k) + g(:,k) - (X_k + G_k) * gamma_k;
        
        if timings
            residual = norm(g(:,k));
        else
            residual =  norm(b-C*x(:,end)); %;
        end
        resvec = [resvec, residual];
        if residual < rel_tol*normb || isnan(residual)
            break
        end
        
        [fx, local_iter] = PMHSS_step(lhs,A,B,b,alpha,V, x(:,k + 1), M0, M1, inner,[], max_inner_iter, x(:,k + 1)); %, residual/2^(k_max-k));
    
        inner_iter = [inner_iter,local_iter];
        g(:,k + 1) = fx - x(:,k + 1);
        % Update increment matrices with new elements.
        X_k = [X_k, x(:,k + 1) - x(:,k)];
        G_k = [G_k, g(:,k + 1) - g(:,k)];
    
        n = size(X_k, 2);
        if n > m_k
            X_k = X_k(:, n - m_k + 1:end);
            G_k = G_k(:, n - m_k + 1:end);
        end
        k = k + 1;
        
    end
    k = k - 1;
    relresidual = residual/normb; 
    xout = x(:,end);
end 


function [xk, ITER] = PMHSS_step(lhs,A,B,b,alpha,V,xk,M0,M1,inner,tol,max_inner_iter,fx)
    if ~exist('tol','var')
         tol = 1e-12;
    end
     
    if all([alpha ~= 1])
        A1 = (alpha*V+A);
        b1 = (alpha*V-1i*B)*xk+b;
        
        [xk05, ~, ITER2] = pcg(A1,b1,1e-13,max_inner_iter,[],[],xk);
        
        A2 = (alpha*V+B);
        b2 = (alpha*V+1i*A)*xk05-1i*b;
        
        [xk, ~, ~, ITER1] = pcg(A2,b2,1e-13,max_inner_iter,[],[],xk05);
        iter_inner = ITER1 + ITER2;
    else
    switch inner 
        case {'\', 'backslash', 'BACKSLSASH'}
            iter_inner = NaN;
            xk = lhs\((1+1i)*(A-1i*B)*xk + (1-1i)*b);
        case {'CG', 'cg', 'PCG', 'pcg'}
            rhs = ((A-1i*B)*((1+1i)*xk) + (1-1i)*b);
            [xk, ~, ~, iter_inner, resvec] = pcg(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk);
        case {'BCG', 'bcg', 'BICG','bicg'}
            rhs = ((1+1i)*(A-1i*B)*xk + (1-1i)*b);
            [xk, ~, ~, iter_inner, resvec] = bicg(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk);
        case {'MINRES', 'minres'}
            rhs = ((1+1i)*(A-1i*B)*xk + (1-1i)*b);
            [xk, ~, ~, iter_inner, ~] = minres(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk);
        case {'SYMMLQ', 'symmlq'}
            rhs = (1+1i)*((A*xk)-(1i*B*xk) -1i*b);
            [xk, ~, ~, iter_inner, ~] = symmlq(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk);
        case {'QMR', 'qmr'}
            rhs = (1+1i)*((A*xk)-(1i*B*xk) -1i*b);
            [xk, ~, ~, iter_inner, resvec] = qmr(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk);
        case {'GMRES', 'gmres'}
            rhs = (1+1i)*((A*xk)-(1i*B*xk) -1i*b);
            [xk, ~, ~, iter_inner, ~] = gmres(lhs,rhs,[],1e-12,max_inner_iter,M0,M1,xk);
        otherwise
            iter_inner = NaN;
            xk = lhs\((1+1i)*(A-1i*B)*xk + (1-1i)*b);
    end
end

ITER = iter_inner;

end
