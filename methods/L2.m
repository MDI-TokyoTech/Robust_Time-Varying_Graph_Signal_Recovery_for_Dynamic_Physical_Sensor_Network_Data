function [Y,S] = L2(L,X_noisy,X,mask,param)
%% set of parameters



epsilon = param.epsilon;
eta = param.eta;
lambda = 1;
maxIter = param.max_iter;
stopping_cri = param.stopping_cri;
phi = @(z) mask'.*z;
phiT = phi;
x = X_noisy';
D = make_D(size(X_noisy,2));


Y = zeros(size(x));
S = zeros(size(x));

y1 = zeros(size(D'*x));
y2 = zeros(size(x));

gamma1 = 0.1;
gamma2 = 0.99/(6*gamma1);

converge_rate = [];
converge_rate1 = [];




for i = 1:maxIter

    Y_pre = Y;
    S_pre = S;
    
    Y = Y - gamma1*(D*y1 + phiT(y2));
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