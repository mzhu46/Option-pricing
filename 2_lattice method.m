%% Jarrow-Rudd Risk Neutral





%Consider pricing a European Call option with the following parameters, 
%X = $60, S0 = $50, r = 5%, ? = 0.2, ?t = 0.01, N = 100.

S0 = 50; X=60; r = 0.05; sigma = 0.2; deltaT = 0.01; N=100;Option_Type='Call';

%Solve the JRR we get the p , d ,u
d = exp((r-sigma^2/2)*deltaT+sigma*deltaT^0.5);
u = exp((r-sigma^2/2)*deltaT-sigma*deltaT^0.5);
E = exp(r*deltaT);

p = (E-d)/(u-d);
%verify the result
%disp(p*u +(1-p)*d - E);

La2 = zeros(N+1);
La2(1,1)=S0; La2(1,2) = S0*u;La2(2,2) = S0*d;

%set upper boundary
for i = 2:N+1
    La2(1,i) = S0*u^(i-1);
end
%generate the lattice
for i = 3:N+1
    for j = 2:i
        La2(j,i) = La2(j-1,i-1)*d;
    end
end

ValueTr = zeros(N+1);

switch Option_Type
    case 'Call'
        for i = 1:N+1
            ValueTr(i,N+1) = max(0,La2(i,N+1)-X);
        end
        
    case 'Put'
        for i = 1:N+1
            ValueTr(i,N+1) = max(0,-La2(i,N+1)+X);
        end
end

for i = N:-1:1 
    for j = 1:N
        ValueTr(j , i) = (p*ValueTr(j,i+1) + (1-p)*ValueTr(j+1,i+1))/E;
    end
end

fprintf("%f4, is call price",ValueTr(1,1))

%% evaluate the JRR and CRR
%X = $60, S0 = $50, r = 5%, ? = 0.2, ?t = 0.01, N = 100.

%The Black-Scholes price for this option is $1.624.

S0 = 50; X=60; r = 0.05; sigma = 0.2; deltaT = 0.01; N=1/deltaT;Option_Type='Call';

a1=EuroOption(S0,X,r,sigma,deltaT,N,Option_Type);
a2=JRRpricing(S0,X,r,sigma,deltaT,N,Option_Type);


for i = 1:4
    deltaT=10^(-i)
    N=floor(1/deltaT)
    a1=EuroOption(S0,X,r,sigma,deltaT,N,Option_Type);
    a2=JRRpricing(S0,X,r,sigma,deltaT,N,Option_Type);
    
    fprintf("CRR method with deltaT = %f   , yields %s price at %f4, ",deltaT,Option_Type,a1-1.624)
    fprintf("\n")
    fprintf("JRR method with deltaT = %f   , yields %s price at %f4, ",deltaT,Option_Type,a2-1.624)
    
end


%JRR method provides closer answer to BS formula when deltaT is large and 
%both of them converges to the theoritical value as deltaT goes smaller

%% Drift comes in

S0 = 50; X=60; r = 0.05; sigma = 0.2; deltaT = 0.005; N=200;Option_Type='Call';

%Solve the JRR we get the p , d ,u, drift ning

sigma2 = 0.2

ning = (log(X)-log(S0))/(deltaT*N)

u = exp(ning*deltaT+sigma*deltaT^0.5);
d = exp(ning*deltaT-sigma2*deltaT^0.5);
E = exp(r*deltaT);
p = (E-d)/(u-d);

La2 = zeros(N+1);
La2(1,1)=S0; La2(1,2) = S0*u;La2(2,2) = S0*d;

%set upper boundary
for i = 2:N+1
    La2(1,i) = S0*u^(i-1);
end
%generate the lattice
for i = 3:N+1
    for j = 2:i
        La2(j,i) = La2(j-1,i-1)*d;
    end
end

ValueTr = zeros(N+1);

switch Option_Type
    case 'Call'
        for i = 1:N+1
            ValueTr(i,N+1) = max(0,La2(i,N+1)-X);
        end
        
    case 'Put'
        for i = 1:N+1
            ValueTr(i,N+1) = max(0,-La2(i,N+1)+X);
        end
end

for i = N:-1:1 
    for j = 1:i
        ValueTr(j , i) = (p*ValueTr(j,i+1) + (1-p)*ValueTr(j+1,i+1))/E;
    end
end

%result = ValueTr(1,1)

%fprintf("%f4, is %s price",ValueTr(1,1),Option_Type)

result = ValueTr(1,1)









