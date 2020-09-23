%This one uses only one loop 


function result = EuroOption2(S0,X,r,sigma,deltaT,N,Option_Type)
%Solve the CRR we get the p , d ,u
d = exp(-sigma*deltaT^0.5);
u = 1/d;
E = exp(r*deltaT);
p = (E-d)/(u-d);
La = zeros(N+1);
La(1,1)=S0; La(1,2) = S0*u;La(2,2) = S0*d;

%set upper boundary
for i = 2:N+1
    La(1,i) = S0*u^(i-1);
end
%generate the lattice
for i = 3:N+1
    La(2:i,i) = La(1:i-1,i-1)*d;
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
        ValueTr(1:N , i) = (p*ValueTr(1:N,i+1) + (1-p)*ValueTr(2:N+1,i+1))/E;
end


result = ValueTr(1,1)
end
