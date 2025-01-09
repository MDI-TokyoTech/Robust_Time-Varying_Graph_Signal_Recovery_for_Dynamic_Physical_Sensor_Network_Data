function [Y,S] = D_X(L,X_noisy,X,mask,param)

%% set of parameters


epsilon = param.epsilon;
eta = param.eta;
maxIter = param.max_iter;
stopping_cri = param.stopping_cri;



phi = @(z) mask'.*z;
phiT = phi;



x = X_noisy';

Y = zeros(size(x));
S = zeros(size(x));
y2 = zeros(size(x));


Lip_L = 10;
Lip = 2*Lip_L;

gamma1 = 2/Lip;
gamma2 = 0.99/ (gamma1*(2));

converge_rate = [];
converge_rate1 = [];
L2 = 2*L;
for i = 1:maxIter

   
    Y_pre = Y;
    S_pre = S;
    
    Y = Y - gamma1*(trace_sum_nabla_2(Y',L2)' + phiT(y2));
    S = ProjFastL1Ball(S - gamma1*(phiT(y2)), eta);
    

    y2_tmp = y2 + gamma2*phi((2*Y-Y_pre)+(2*S-S_pre));
    y2 = y2_tmp - gamma2*ProjFroball(y2_tmp/gamma2, X_noisy',epsilon);
    converge_rate = cat(2, converge_rate, sqrt(sum((Y - Y_pre).^2, 'all')/sum(Y_pre.^2, 'all')));
    converge_rate1 = cat(2, converge_rate1, sqrt(sum((S - S_pre).^2, 'all')/sum(S_pre.^2, 'all')));

    if i>=2 && converge_rate(i) < 10^(-stopping_cri) && converge_rate1(i) < 10^(-stopping_cri)

    end

end

Y = Y';
S = S';
