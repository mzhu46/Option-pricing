% MC method is suitable for Asian Option
%which is path dependent and Asian Option is on the average price 
% instead of the spot price at exercise day

%%INPUYS
S0 =50;       % Price of underlying today
X = 50;       % Strike at expiry
mu = 0.04;    % expected return
sigma = 0.1;    % expected vol.
r = 0.03;     % Risk free rate

dt = 1/365;   % time steps
T = dt*50; % years to expiry
simN = 1000;
E = exp(r*T);


drift = (mu - sigma^2/2)*dt;
sto = sigma*dt^0.5*normrnd(0,1);

PricePath = zeros(simN,round(T/dt)+1);

%initial the PricePath with S0
PricePath(:,1) = S0;
for i = 2:round(T/dt)+1
    for j = 1:simN
        PricePath(j,i) = PricePath(j,i-1) *exp(drift+sigma*dt^0.5*normrnd(0,1));
    end
end

CallPriceM = max(mean(PricePath,2)-X,0);
CallPrice = mean(CallPriceM)/E


%% plot the result
plot(PricePath(1,:))
hold on 
for i = 2:simN
    plot(PricePath(i,:))
end
hold off



%% Matrix method
%%INPUYS
S0 =50;       % Price of underlying today
X = 50;       % Strike at expiry
mu = 0.04;    % expected return
sigma = 0.1;    % expected vol.
r = 0.03;     % Risk free rate
dt = 1/365;   % time steps
T = dt*50; % years to expiry
simN = 10000;
E = exp(r*T);


drift = (mu - sigma^2/2)*dt;
Driftvector = zeros(simN,1);


Driftvector(:,1)=drift;
sto = sigma*dt^0.5*normrnd(0,1,[simN,1]);



PricePath = zeros(simN,round(T/dt)+1);

%initial the PricePath with S0
PricePath(:,1) = S0;
for i = 2: round(T/dt)+1
    %need new sto term each day
    sto = sigma*dt^0.5*normrnd(0,1,[simN,1]);
    %.* one by one
    PricePath(:,i) = PricePath(:,i-1).*exp(Driftvector+sto);  
end
CallPriceM = max(mean(PricePath,2)-X,0);
CallPriceM2 = max(PricePath(:,round(T/dt)+1)-X,0);
OptionPrice2European = mean(CallPriceM2)/E
CallPriceAsian = mean(CallPriceM)/E
