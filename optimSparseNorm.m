function f = optimSparseNorm(g,T,lambda,k)
% optimizes ||g - Tf||_2^2 + lambda^2 || D|f| ||_k^k
fo = minres(T'*T,T'*g);





f = fo;
end

function out = Hessian(f,T,lambda2,D,k)
    out = 2*T'*T + k*lambda2^2*Phi(f)'*...
        D'*Lambda2(f,D,k,eps)*D*Phi(f);
end

function out = Lambda2(f,D,k,eps)
    grads = D*abs(f);
    weights = 1./(grads.^2 + eps).^(1-k/2);
    out = diag(weights);
end

function out = Phi(f)
    out = diag(exp(-1j*angle(f)));
end