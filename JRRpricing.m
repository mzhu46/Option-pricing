function JRRprice = JRR(S0,X,r,sigma,deltaT,N,Option_Type)

d = exp((r-sigma^2/2)*deltaT+sigma*deltaT^0.5);
u = exp((r-sigma^2/2)*deltaT-sigma*deltaT^0.5);
E = exp(r*deltaT);

p = (E-d)/(u-d);
%verify the result
%disp(p*u +(1-p)*d - E);

La = zeros(N+1);
La(1,1)=S0; La(1,2) = S0*u;La(2,2) = S0*d;

%set upper boundary
for i = 2:N+1
    La(1,i) = S0*u^(i-1);
end
%generate the lattice
for i = 3:N+1
    for j = 2:i
        La(j,i) = La(j-1,i-1)*d;
    end
end

ValueTr = zeros(N+1);

switch Option_Type
    case 'Call'
        for i = 1:N+1
            ValueTr(i,N+1) = max(0,La(i,N+1)-X);
        end
        
    case 'Put'
        for i = 1:N+1
            ValueTr(i,N+1) = max(0,-La(i,N+1)+X);
        end
end

for i = N:-1:1 
    for j = 1:N
        ValueTr(j , i) = (p*ValueTr(j,i+1) + (1-p)*ValueTr(j+1,i+1))/E;
    end
end
%fprintf("%f4, is %s price",ValueTr(1,1),Option_Type)

JRRprice = ValueTr(1,1)
end