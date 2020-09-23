% Script to price an Asian put option using a monte-carlo approach.
S0 = [50 48]  ;       % Price of underlying today
mu = [0.03 0.06];     % expected return
sig = [0.05 0.1];     % expected vol.
corr = [1 0.25;0.25 1]; % correlation matrix
dt = 1/100;   % time steps
etime = 100;   % days to expiry
T = dt*etime; % years to expiry

nruns = 1000; % Number of simulated paths

% Generate potential future asset paths
S = AssetPathsCorrelated(S0,mu,sig,corr,dt,etime,nruns);

%% Check the correlation
simN = nruns
R1 = zeros(simN,2,2);
for K1 = 1:simN

    A = [S(:,K1,1) , S(:,K1,2)];

    R1(K1,:,:) = corrcoef((A));
end

mean(R1(:,1,2))
mean(R1(:,2,1))

mean(R1,1)

%% check 
result = zeros(2,1);
for p = 0.01:0.01:0.9
corr = [1 p;p 1]; % correlation matrix
S = AssetPathsCorrelated(S0,mu,sig,corr,dt,etime,nruns);

simN = nruns;
R1 = zeros(simN,2,2);
for K1 = 1:simN

    A = [S(:,K1,1) , S(:,K1,2)];

    R1(K1,:,:) = corrcoef((A));
end

mean(R1(:,1,2));
diff = [ p ; mean(R1(:,1,2))];

result = [result diff];
end
%% plot the result
plot(result(1,:))
hold on 
plot(result(2,:))
hold off 


function S = AssetPathsCorrelated(S0,mu,sig,corr,dt,steps,nsims)
% Function to generate correlated sample paths for assets assuming
% geometric Brownian motion.
%
% S = AssetPathsCorrelated(S0,mu,sig,corr,dt,steps,nsims)
%
% Inputs: S0 - stock price
%       : mu - expected return
%       : sig - volatility
%       : corr - correlation matrix
%       : dt - size of time steps
%       : steps - number of time steps to calculate
%       : nsims - number of simulation paths to generate
%
% Output: S - a (steps+1)-by-nsims-by-nassets 3-dimensional matrix where
%             each row represents a time step, each column represents a
%             seperate simulation run and each 3rd dimension represents a
%             different asset.
%
% Notes: This code focuses on details of the implementation of the
%        Monte-Carlo algorithm.
%        It does not contain any programatic essentials such as error
%        checking.
%        It does not allow for optional/default input arguments.
%        It is not optimized for memory efficiency or speed.

% Author: Phil Goddard (phil@goddardconsulting.ca)
% Date: Q2, 2006

% get the number of assets
nAssets = length(S0);

% calculate the drift
nu = mu - sig.*sig/2;

% do a Cholesky factorization on the correlation matrix
R = chol(corr);
% pre-allocate the output
S = nan(steps+1,nsims,nAssets);

% generate correlated random sequences and paths
for idx = 1:nsims
    % generate uncorrelated random sequence
    x = randn(steps,size(corr,2));
    % correlate the sequences
    ep = x*R;

    % Generate potential paths
    S(:,idx,:) = [ones(1,nAssets); ...
        cumprod(exp(repmat(nu*dt,steps,1)+ep*diag(sig)*sqrt(dt)))]*diag(S0);
end

% If only one simulation then remove the unitary dimension
if nsims==1
    S = squeeze(S);
end   
end