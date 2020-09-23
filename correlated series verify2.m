S0 = [50 48]  ;       % Price of underlying today
mu = [0.03 0.06];     % expected return
sig = [0.05 0.1];     % expected vol.
corr = [1 0.25;0.25 1]; % correlation matrix
dt = 1/365;   % time steps
etime = 50;   % days to expiry
T = dt*etime; % years to expiry

nruns = 10000; % Number of simulated paths

mu1 = 0.03 - 1/2 * 0.05^2 
sigma1 = 0.05
%%  try

path = GBMsimulation(50,0.03,0.05,1/365,50,1000);
















function GBM_path = GBMsimulation ( S0, mu ,sigma ,dt, steps,simN)
%calculate the drift term
drift = (mu - sigma^2/2)*dt;
Driftvector = zeros(simN,1);
Driftvector(:,1)=drift;
sto = sigma*dt^0.5*normrnd(0,1,[simN,1]);


PricePath = zeros(simN,round(steps)+1);
%initial the PricePath with S0
PricePath(:,1) = S0;
%Get Price Paths
for i = 2: round(steps)+1
    %need new sto term each day
    sto = sigma*dt^0.5*normrnd(0,1,[simN,1]);
    %.* one by one
    PricePath(:,i) = PricePath(:,i-1).*exp(Driftvector+sto);  
end
GBM_path =PricePath;

end
