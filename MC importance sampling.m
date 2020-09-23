%% counting example 
format long
for i = 1:100
   a(i) = sum(randn(1000,1)>=3)/100;
end
meanMC  = mean(a)
varMC   = var(a)

%% IS
for j = 1:100
   N = 100;
   x = randn (N,1) + 4;
     for ii = 1:N
         h = x(ii)>=3;
         b = exp(8-4*x(ii));
         w(ii) = h*b;
     end
  I(j) = sum(w)/N;
end
MEAN = mean(I)
VAR = var(I)
%% IS with pdf

for j = 1:100
   N = 100;
   x = randn (N,1) + 4;
     for ii = 1:N
         h = x(ii)>=3;
         b = normpdf(x(ii),0,1)/normpdf(x(ii),4,1);
         w(ii) = h*b;
     end
  I(j) = sum(w)/N;
end
MEAN = mean(I)
VAR = var(I)



%% pricing example
S0 =50;       % Price of underlying today
X = 80;       % Strike at expiry
sigma = 0.4;    % expected vol.
r = 0.05;     % Risk free rate
dt = 1/1000;   % time steps
T = 1; % years to expiry
simN = 100000;
E = exp(r*T);


N=T/dt;
% deep out of money option price matrix
dt = dt*N;
alpha = (r- sigma^2/2)*dt;
beta = (log(X/S0)/T -sigma^2/2)*dt;
sig=sigma*dt^0.5;
Z_Random = normrnd(0,1,[simN,1]);

A1 = S0*exp(alpha + sig*Z_Random);

hist(A1)
A1_payoff = zeros(simN,1);
for i = 1: simN
    A1_payoff(i) = max(A1(i)-X,0);
end 
[Price, VarPrice, CI] = normfit(A1_payoff)

%%%%IS
A2 = S0*exp(beta + sig*Z_Random);

% histogram(A2,50)
% hold on 
% histogram(A1,50)
% hold off

% get h(x) and ratio = f(x)/g(x) and f and g are pdf of norm

A2_payoff = zeros(simN,1);
for i = 1: simN
    A2_payoff(i) = max(A2(i)-X,0);
end 
%% check the distribution
Es = S0*exp((beta+sigma^2/2)*T)
mu = log(X/S0)/T;
Es2 = S0*exp(mu*T)
Evar= S0^2*exp(2*mu*T)*(exp(sigma^2*T)-1)
var(A2)
%[Price, VarPrice, CI] = normfit(A2)
%ratio 


%%%%%
%%
mu1 = 50 ; sigma1 = 20;
mu2 = 80 ; sigma2 = 25;
N=2000
X1 = 50+ sigma1 * randn(N,1);
K=80
for j = 1:100
   N = 100;
   x = mu2 + sigma2*randn(N,1);
     for ii = 1:N
         h = x(ii)>=K;
         b = normpdf(x(ii),mu1,sigma1)/normpdf(x(ii),mu2,sigma2);
         B(ii)=b;
         w(ii) = h*b;
     end
  I(j) = sum(w)/N;
end
MEAN = mean(I)
mean2 = mean(X1>K)
VAR = var(I)
VAR2 = var(X1>K)
%% apply on pricing

mu_f = mean(A1); var_f = var(A1)^0.5;
mu_g = mean(A2); var_g = var(A2)^0.5;


K=80
for j = 1:100
   N = 1000;
   x = mu_g + var_g*randn(N,1);
     for ii = 1:N
         h = max(x(ii)-K,0);
         b = normpdf(x(ii),mu_f,var_f)/normpdf(x(ii),mu_g,var_g);
         B(ii)=b;
         w(ii) = h*b;
     end
  I(j) = sum(w)/N;
end

mean(I)
mean(max(A1-K,0))