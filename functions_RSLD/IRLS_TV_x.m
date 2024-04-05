
% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(ux)
% where ux is the first of three images in u, which has size MxNx3
function u = IRLS_TV_x(b,A,mu,M,N,tol,minimask)

AtA = A'*A;
Atb = A'*b;

%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dy = kron(speye(N),D);
Dy = kron(sparse(1,1,1,3,3), Dy);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dx = kron(D,speye(M));
Dx = kron(sparse(1,1,1,3,3), Dx);

D = [Dx' Dy']';
% 
ite_irls = 0;
error = 1;

%[u,~] = cgs(AtA+mu*(D')*D,Atb);
[u,~] = cgs(AtA,Atb);
G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u(1:M*N),M,N,minimask(1:M*N));

while error > tol && ite_irls < 200
    
    ite_irls = ite_irls + 1;
    Dh = Dx*u;
    Dv = Dy*u;
    vksquare = Dh.^2 + Dv.^2;
    
    eps = 0.3;
    P = sqrt(vksquare + eps^2);
    P = 1./P;
    %Huber seems not to work;
    
    %P = sqrt(P.^2 + eps^2);
    P = P(:).*minimask;
    omega = speye(M*N*3);   % sparse diagonal of just ones
    omega = spdiags(P,0,omega);   % sparse diagonal of P values instead of ones.
    W = kron(speye(2),omega);
    
    [u,~] = cgs(AtA + mu*D'*W*D, Atb, 1e-6, 200);
    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u(1:M*N),M,N,minimask(1:M*N));
    error = abs(G(ite_irls+1) - G(ite_irls));
    
end

%figure(909); plot(1:length(G),G);

end
