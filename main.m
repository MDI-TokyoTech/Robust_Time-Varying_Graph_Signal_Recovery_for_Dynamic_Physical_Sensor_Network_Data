function main()


addpath('dataset')
addpath('methods')
addpath('utils')
load("synthetic_dataset.mat")

rng('default')
profile on;


close all;

%%set parameters
nreplicate = 5;
stopping_cri = 4;
max_iter = 5000;
sigma = 0.1;
p_s = 0.1;
p_p = 0.1;
lambda = 1;


N = size(X,1);
M = size(X,2);
epsilon = 0.9*sqrt(N*M*(1-p_s)*(1-p_p))*sigma;
eta = 0.9*p_s*0.5*N*M;

param = struct();
param.stopping_cri = stopping_cri;
param.max_iter = max_iter;
param.sigma = sigma;
param.p_s = p_s;
param.p_p = p_p;
param.N = N;
param.M = M;
param.epsilon = epsilon;
param.eta = eta;
param.lambda = lambda;
num_of_methods = 5;

result_rmse = zeros(num_of_methods,nreplicate);
result_mae = zeros(num_of_methods,nreplicate);

original_rmse = zeros(nreplicate,1);
original_mae = zeros(nreplicate,1);



signal = zeros(N,M, nreplicate);
signal_noisy = zeros(N,M,nreplicate);

laplacian = zeros(N,N,M,nreplicate);
graphs = cell(nreplicate,1);
Y = zeros(N,M,num_of_methods,nreplicate);
S_bar = zeros(N,M,nreplicate);

dir_name = ['results/',num2str(N),'_',num2str(M)];
name = [dir_name,'/sigma=',num2str(sigma),'_P_s=',num2str(p_s),'_P_p=',num2str(p_p),'_result.mat'];



for ii = 1:nreplicate

tmp = rand(N,M);
tmp = tmp<p_s;
S_bar(:,:,ii) = (-1 + 2*rand(N,M)).*tmp;
mask = rand(N,M);
mask = mask>p_p;
phi = @(z) mask.*z;

X_noisy = phi(X + S_bar(:,:,ii) + sigma*randn(size(X)));

graphs{ii} = graph;

signal(:,:,ii) = X;
signal_noisy(:,:,ii) = X_noisy;

laplacian(:,:,:,ii) = L;


method = 1;
[Y(:,:,method,ii),~] = L2(L(:,:,M/2),X_noisy,X,mask,param);
result_rmse(method,ii) = rms(X-Y(:,:,method,ii),'all');
result_mae(method,ii) = mean(abs(X-Y(:,:,method,ii)),'all');

method = 2;
[Y(:,:,method,ii),~] = S_X(L(:,:,M/2),X_noisy,X,mask,param);
result_rmse(method,ii) = rms(X-Y(:,:,method,ii),'all');
result_mae(method,ii) = mean(abs(X-Y(:,:,method,ii)),'all');

method = 3;
result_rmse(method,ii) = rms(X-Y(:,:,method,ii),'all');
result_mae(method,ii) = mean(abs(X-Y(:,:,method,ii)),'all');

method = 4;
[Y(:,:,method,ii),~] = S_X_L2(L(:,:,M/2),X_noisy,X,mask,param);
result_rmse(method,ii) = rms(X-Y(:,:,method,ii),'all');
result_mae(method,ii) = mean(abs(X-Y(:,:,method,ii)),'all');

method = 5;
[Y(:,:,method,ii),~] = D_X_L2(L,X_noisy,X,mask,param);
result_rmse(method,ii) = rms(X-Y(:,:,method,ii),'all');
result_mae(method,ii) = mean(abs(X-Y(:,:,method,ii)),'all');


original_rmse(ii) = rms(X-X_noisy,'all');
original_mae(ii) = mean(abs(X-X_noisy),'all');

end



if exist(dir_name,'dir') == 0
    mkdir (dir_name)
end
save(name,'-v7.3');

