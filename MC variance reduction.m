% Antithetic Sampling is generating a negative correlated path for each
% existing path to do variance reduction by negative cov

%set parameter
%(S0,X,mu,r,sigma,dt,T,Option_Type,simN)
%(50,50,0.04,0.03,0.1,1/365,50/365,'Call',100)

%sigma bar = 0.6735
sigma_par=MCsimulation_AsianOption_Matrix(50,50,0.04,0.03,0.1,1/365,50/365,'Call',100000);
var(sigma_par)^0.5

%
%We can easily verify that BY the CLT the sample variance get reduced by
%dividing N^0.5
p2 = zeros(100,1);
N=30;
for i = 1:1000
    p2(i,1) = MCsimulation_AsianOption(50,50,0.04,0.03,0.1,1/365,50/365,'Call',N);
end
var(p2)^0.5  * N^0.5







%% EuroOption

%sigma bar = 1.1985
sigma_par=MCsimulation_EuroOption_Matrix(50,50,0.04,0.03,0.1,1/365,50/365,'Call',100000);
var(sigma_par)^0.5;

[Price, VarPrice, CI]=normfit(sigma_par)


%% draw negative correlated paths

%%INPUYS
S0 =50;       % Price of underlying today
X = 50;       % Strike at expiry
mu = 0.04;    % expected return
sigma = 0.1;    % expected vol.
r = 0.03;     % Risk free rate
dt = 1/365;   % time steps
T = dt*50; % years to expiry
simN = 100000;
E = exp(r*T);


drift = (mu - sigma^2/2)*dt;
Driftvector = zeros(simN,1);


Driftvector(:,1)=drift;
sto = sigma*dt^0.5*normrnd(0,1,[simN,1]);



PricePath = zeros(simN,round(T/dt)+1);

%initial the PricePath with S0
PricePath(:,1) = S0;
PricePath2= PricePath;
for i = 2: round(T/dt)+1
    %need new sto term each day
    %now only generate half random
    sto = sigma*dt^0.5*normrnd(0,1,[simN,1]);
    %.* one by one
    PricePath(:,i) = PricePath(:,i-1).*exp(Driftvector+sto); 
    PricePath2(:,i) = PricePath2(:,i-1).*exp(Driftvector-sto); 
end
CallPriceM = max(PricePath(:,round(T/dt)+1)-X,0);
CallPriceM2 = max(PricePath2(:,round(T/dt)+1)-X,0);
%doing averaging 
CallPriceM3 = 0.5*(CallPriceM+CallPriceM2);
[Price, VarPrice, CI] = normfit(CallPriceM3)

%% verify the variance formula
Cov(CallPriceM,CallPriceM2)
var(CallPriceM3)


%%%%
%% Important Sampling
% Initalize Input Variables
randn('state',0)
Option_Type='Call';
S0 =50;       % Price of underlying today
X = 80;       % Strike at expiry
sigma = 0.4;    % expected vol.
r = 0.05;     % Risk free rate
dt = 1/1000;   % time steps
T = 1; % years to expiry
simN = 100000;
E = exp(r*T);

N=T/dt

% deep out of money option price matrix

dt = dt*N
alpha = (r- sigma^2/2)*dt;
beta = (log(X/S0)/T -sigma^2/2)*dt;
sig=sigma*dt^0.5;
Z_Random = normrnd(0,1,[simN,1]);
VY = beta + sig*Z_Random;
VX = alpha + sig*Z_Random;
Weights = exp( (beta^2 - alpha^2 + 2*(alpha - beta)*VY)/(2*sig^2));



%%%%%%%
E = exp(r*T);

Driftvector = zeros(simN,1);
Driftvector2 =Driftvector; 
Driftvector(:,1)=alpha;
Driftvector2(:,1)=beta;


PricePath = zeros(simN,round(T/dt)+1);
%initial the PricePath with S0
PricePath(:,1) = S0;

PricePath2 = PricePath;
%Get Price Paths
for i = 2: round(T/dt)+1
    %need new sto term each day
    Z_Random = normrnd(0,1,[simN,1]);
    VY = beta + sig*Z_Random;
    sto = sigma*dt^0.5*Z_Random;
    %generate weights
    weight=exp( (beta^2 - alpha^2 + 2*(alpha - beta)*VY)/(2*sig^2));
    %.* one by one
    PricePath(:,i) = PricePath(:,i-1).*exp(Driftvector+sto);  
    PricePath2(:,i) = PricePath2(:,i-1).*exp((Driftvector2+sto));  
end

switch Option_Type
    case 'Call'
        CallPriceM = max(PricePath(:,round(T/dt)+1)-X,0);
        CallPriceM2 = max(PricePath2(:,round(T/dt)+1)-X,0).*weight;
        OptionPrice = mean(CallPriceM)/E;
    case 'Put'
        PutPriceM = max(X-mean(PricePath,2),0);
        OptionPrice = mean(PutPriceM)/E;
end
%fprintf("%f4, is %s price",OptionPrice,Option_Type)
result = CallPriceM;
result_IS = CallPriceM2;


[Price, VarPrice, CI] = normfit(result)
[PriceIS, VarPriceIS, CIIS] = normfit(result_IS)
