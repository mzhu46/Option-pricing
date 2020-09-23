for i = 1:1000000
   a(i) = sum(randn(100,1)>=3)/100;
end
meanMC  = mean(a)
varMC   = var(a)

%% importance sampling
Var_V = zeros(10,1);
for i = 2:12
    for j = 1:1000
       N = 1000;
       x = randn (N,1) + i;
         for ii = 1:N
             h = x(ii)>=3;
             b = exp(i^2/2-i*x(ii));
             w(ii) = h*b;
         end
      I(j) = sum(w)/N;
    end
    Var_V(i,1) = var(I);
end
plot(Var_V)



%
i=10;
for j = 1:1000
   N = 1000;
   x = randn (N,1) + i;
     for ii = 1:N
         h = x(ii)>=3;
         b = exp(i^2/2-i*x(ii));
         w(ii) = h*b;
     end
  I(j) = sum(w)/N;
end
var(I)



%%
VAR_V = zeros(50,2);
for i = 2:0.2:6
    for j = 1:1000
       N = 1000;
       x = randn (N,1) + i;
         for ii = 1:N
             h = x(ii)>=3;
             b = exp(i^2/2-i*x(ii));
             w(ii) = h*b;
         end
      I(j) = sum(w)/N;
    end
    VAR_V(round((i-2)*5+1),1) = mean(I);
    VAR_V(round((i-2)*5+1),2) = var(I);
end
plot(VAR_V(:,2))

%% study the drift of Importance Sampling
%% Importance Sampling (Out of Money)

clear;
clc;
randn('state',0)

% Initalize Input Variables
S0 = 50;
K = 80;
r = 0.05;
T = 5/12;
sigma = 0.4;
NRepl = 1000000;


Var_changewithdrift = zeros(30,1);

for i = 1:100
    
    % Calculate Model Parameters
    nuT = (r - 0.5*sigma^2)*T;
    siT = sigma * sqrt(T);
    ISnuT = log((K-10+i)/S0) - 0.5*sigma^2*T;
    Veps = randn(NRepl, 1); % generating the noise for the dz
    VY = ISnuT + siT*Veps;



    ISRatios = exp((2 * (nuT - ISnuT) * VY - nuT^2 + ISnuT^2)/2/siT^2);

    % Calculating payoffs
    Payoff = max(0, S0*exp(nuT+siT*Veps) - K);
    PayoffIS = max(0, (S0*exp(VY)-K)) .* ISRatios; % notice the adustment for IS

    DiscPayoff = exp(-r*T) * (Payoff);
    DiscPayoffIS = exp(-r*T)*PayoffIS;

    [Price, VarPrice, CI] = normfit(DiscPayoff);
    [PriceIS, Var_changewithdrift(i,1), CIIS] = normfit(DiscPayoffIS);
end

plot(Var_changewithdrift)