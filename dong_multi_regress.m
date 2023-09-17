function [beta,ResMS,T,P] = dong_multi_regress(Y,X,TR)
% Multiple regression
% Input:
%      Y: Dependent variable: n1 samples x N dimensions
%      X: Independent variable: n1 samples x M dimensions
%      TR: fMRI repeating time (optional).If input TR, it will removes low 
%         frequency confounds X0 for fMRI data.
% Output:
%      beta: regression coeffeients. M dimensions X N dimensions
%      T: t-values of regression coeffeients. M dimensions X N dimensions
%      P: p-values of regression coeefients. M dimensions X N dimensions
%      ResMS: the variance of the error.   1 X N dimensions
% -------------------------------------------------------------------------
% Code Summary for working in School of Life Science and Technology,UESTC.
% This software is for non commercial use only.
% It is freeware but not in the public domain.
% -------------------------------------------------------------------------
% copyright is Li Dong
% Email: Li_dong729@163.com
% 12.15.2014 last edit.
% Revised by Li Dong (2015.5.1): 1) change SSE calculating method to save
%                                memory.
%                                2) using pinv to calculate inverse operation
% -------------------------------------------------------------------------
if nargin < 2
    error('Reqiured at least 2 inputs !');
elseif nargin == 2
    TR = [];
end

N = size(Y,2);
M = size(X,2);
T = zeros(M,N);
P = zeros(M,N);
[nrow,ncol] = size(X);
% -------------------------
% Removes low frequency confounds X0 for fMRI data
if ~isempty(TR)
    K.RT = TR;
    K.row = 1:nrow;
    K.HParam = 128; % cut-off period in seconds. Default is 128s in SPM8
    [K] = spm_filter(K);
    Y = spm_filter(K,Y);
    X = spm_filter(K,X);
end
% -------------------------
beta = X \ Y; % regression coeffeients
% beta_std = beta.*(repmat(std(X)',1,N)./repmat(std(Y),M,1)); % standard regression coeffeients???
r = Y - X * beta;
SSE = sum(r.^2);
ResMS = SSE ./ (nrow-ncol);
if nargout >= 4
    for i = 1:N
        covb = ResMS(i) * pinv(X'*X);
        T_all =  beta(:,i)./sqrt(diag(covb));
        T_all(~isfinite(T_all)) = 0; % Toss NaN
        T(:,i) = T_all(1:M);
        P(:,i) = 2*tcdf((-abs(T(:,i))),(nrow-ncol));
    end
end
