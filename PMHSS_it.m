% function [] = PMHSS_it(C,x0,b,rel_tol,k_max,M0,M1,inner)
function [xout, k, inner_iter, resvec, relresidual] = PMHSS_it(C,xk,b,rel_tol,k_max,M0,M1,inner,max_inner_iter,timings)
    bnorm = norm(b);
    A = real(C);
    B = imag(C);
    lhs = 2*(A+B);
    V = A;
    alpha = 1; 
    inner_iter = []; resvec = [norm(C*xk-b)];
    for k = 1:k_max
        xk2 = xk; 
        [xk, ITER] = PMHSS_step(lhs,A,B,b,alpha,V,xk,M0,M1,inner,1e-12,max_inner_iter);
        if timings 
            residual = norm(xk2-xk);
        else 
            residual = norm(C*xk-b);
        end
        inner_iter = [inner_iter, ITER];
        resvec = [resvec, residual];
        if residual < rel_tol*bnorm
            break
        end
        
    end
    
    xout = xk;
    relresidual = residual/bnorm;

end

function [xk, ITER] = PMHSS_step(lhs,A,B,b,alpha,V,xk,M0,M1,inner,tol,max_inner_iter)
    
    %1-1i; %-1i+0.01; %-1i-0.01; %-1*1i-0.01;
    if alpha ~= 1
        A1 = (alpha*V+A);
        b1 = (alpha*V-1i*B)*xk+b;
        
        [xk05, ~, ITER2] = pcg(A1,b1,1e-13,200,[],[],xk); 
        %xk05 = (A1\b1);
        
        A2 = (alpha*V+B);
        b2 = (alpha*V+1i*A)*xk05-1i*b;
        
        [xk, ~, ~, ITER1] = pcg(A2,b2,1e-13,200,[],[],xk05); %
        %xk = (A2\b2);
        iter_inner = ITER1 + ITER2;
    else
        switch inner 
            case {'\', 'backslash', 'BACKSLSASH'}
                iter_inner = NaN;
                xk = lhs\((1+1i)*(A-1i*B)*xk + (1-1i)*b);
            case {'CG', 'cg', 'PCG', 'pcg'}
                rhs = ((A-1i*B)*((1+1i)*xk) + (1-1i)*b);
                [xk, ~, ~, iter_inner, resvec] = pcg(lhs,rhs,tol,max_inner_iter,M0,M1,xk);
            case {'BCG', 'bcg', 'BICG','bicg'}
                %(A+B),1/2*((1+1i)*(A-1i*B)*xk + (1-1i)*b)
                rhs = ((1+1i)*(A-1i*B)*xk + (1-1i)*b);
                [xk, ~, ~, iter_inner, ~] = bicg(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk);
            case {'MINRES', 'minres'}
                rhs = ((1+1i)*(A-1i*B)*xk + (1-1i)*b);
                [xk, ~, ~, iter_inner, ~] = minres(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk);
            case {'SYMMLQ', 'symmlq'}
                rhs = (1+1i)*((A*xk)-(1i*B*xk) -1i*b);
                [xk, ~, ~, iter_inner, ~] = symmlq(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk);
            case {'QMR', 'qmr'}
                rhs = (1+1i)*((A*xk)-(1i*B*xk) -1i*b);
                [xk, ~, ~, iter_inner, ~] = qmr(lhs,rhs,1e-12,max_inner_iter,M0,M1,xk); 
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