function [Y,S] = D_X_L2(L,X_noisy,X,mask,param)

%% set of parameters

epsilon = param.epsilon;
eta = param.eta;
lambda = param.lambda;
maxIter = param.max_iter;
stopping_cri = param.stopping_cri;
phi = @(z) mask'.*z;
phiT = phi;
x = X_noisy';
D = make_D(size(X_noisy,2));


lambda = lambda*(trace(X'*L(:,:,1)*X) / norm(D'*X',2));

Y = zeros(size(x));
S = zeros(size(x));

y1 = zeros(size(D'*x));
y2 = zeros(size(x));


Lip_L = 10;
Lip = 2*Lip_L;

gamma1 = 1/Lip;
gamma2 = 0.49/ (gamma1*(6));


converge_rate = [];
converge_rate1 = [];




for i = 1:maxIter
   
    Y_pre = Y;
    S_pre = S;
    
    Y = Y - gamma1*(trace_sum_nabla(Y',L)' + D*y1 + phiT(y2));
    S = ProjFastL1Ball(S - gamma1*(phiT(y2)), eta);
    
    y1_tmp = y1 + gamma2*D'*(2*Y-Y_pre);
    y2_tmp = y2 + gamma2*phi((2*Y-Y_pre)+(2*S-S_pre));
    y1 = y1_tmp - gamma2*ProxL2norm(y1_tmp/gamma2, lambda/gamma2);
    y2 = y2_tmp - gamma2*ProjFroball(y2_tmp/gamma2, X_noisy',epsilon);
    
    converge_rate = cat(2, converge_rate, sqrt(sum((Y - Y_pre).^2, 'all')/sum(Y_pre.^2, 'all')));
    converge_rate1 = cat(2, converge_rate1, sqrt(sum((S - S_pre).^2, 'all')/sum(S_pre.^2, 'all')));
    
    if i>=2 && converge_rate(i) < 10^(-stopping_cri) && converge_rate1(i) < 10^(-stopping_cri)
         break
    end
end
Y = Y';
S = S';