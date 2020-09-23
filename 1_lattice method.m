% CRR method    
% Matching the return
%p*u + (1-p)*d = e^(r*deltaT)

%Matching the variance
%Var(X) = E(X^2) - E(X)^2
%sigma^2*deltaT = p*u^2 + (1-p)*d^2 - (e^(r*deltaT))^2

%Cox-Ross-Rubinstein
%we get another equation
%Assume u =1/d to recommbine the trees to get less nodes

%%
%Consider pricing a European Call option with the following parameters, 
%X = $60, S0 = $50, r = 5%, ? = 0.2, ?t = 0.01, N = 100.

S0 = 50; X=60; r = 0.05; sigma = 0.2; deltaT = 0.01; N=100;Option_Type='Call';

%Solve the CRR we get the p , d ,u
d = exp(-sigma*deltaT^0.5);
u = 1/d;
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

fprintf("%f4, is call price",ValueTr(1,1))



%%
%American Option can exercise early



S0 = 50; X=60; r = 0.05; sigma = 0.2; deltaT = 0.01; N=100;Option_Type='Call ';

%Solve the CRR we get the p , d ,u
d = exp(-sigma*deltaT^0.5);
u = 1/d;
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
        
        
        for i = N:-1:1 
            for j = 1:i
                ValueTr(j , i) =max(La(j,i)-X, (p*ValueTr(j,i+1) + (1-p)*ValueTr(j+1,i+1))/E);
            end
        end

        
        
        
    case 'Put'
        for i = 1:N+1
            ValueTr(i,N+1) = max(0,-La(i,N+1)+X);
        end
        
        for i = N:-1:1 
            for j = 1:i
                ValueTr(j , i) =max(X-La(j,i), (p*ValueTr(j,i+1) + (1-p)*ValueTr(j+1,i+1))/E);
            end
        end
end



fprintf("  %f4, is %s price",ValueTr(1,1),Option_Type)




