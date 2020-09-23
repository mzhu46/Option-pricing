function result = MCsimulation_AsianOption_Matrix(S0,X,mu,r,sigma,dt,T,Option_Type,simN)

%get necessary term 
%E discounting factor and drift term
E = exp(r*T);
drift = (mu - sigma^2/2)*dt;
Driftvector = zeros(simN,1);
Driftvector(:,1)=drift;
sto = sigma*dt^0.5*normrnd(0,1,[simN,1]);

PricePath = zeros(simN,round(T/dt)+1);
%initial the PricePath with S0
PricePath(:,1) = S0;
%Get Price Paths
for i = 2: round(T/dt)+1
    %need new sto term each day
    sto = sigma*dt^0.5*normrnd(0,1,[simN,1]);
    %.* one by one
    PricePath(:,i) = PricePath(:,i-1).*exp(Driftvector+sto);  
end

switch Option_Type
    case 'Call'
        CallPriceM = max(mean(PricePath,2)-X,0);
        OptionPrice = mean(CallPriceM)/E;
        
    case 'Put'
        PutPriceM = max(X-mean(PricePath,2),0);
        OptionPrice = mean(PutPriceM)/E;
end 
result = CallPriceM;
end