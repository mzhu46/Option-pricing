K = 60;S0=50;r=0.05;sigma=0.2;Svec=0:0.1:100;tvec=0:0.001:1;
tic
deltat = tvec(2)-tvec(1);
ds = Svec(2)-Svec(1);
N = length(Svec)-1;
M = length(tvec)-1;
% j is always follow the S vector pace
j = 0:N;
sig2 = sigma*sigma;
aj = (deltat*j/2).*(r - sig2*j);
bj = 1 + deltat*(sig2*(j.^2) + r);
cj = -(deltat*j/2).*(r + sig2*j);

Bn = diag(aj(3:N),-1) + diag(bj(2:N)) + diag(cj(2:N-1),1);



sol2=zeros(N+1,M+1);
%set value and set boundary first row to be 0 and last row be maximum value
for i = 1:N+1
    sol2(i,1) = max(i-1-K,0);
end
for j = 1:M+1
    sol2(N+1,j)=sol2(N+1,1);
end
for t = 1:M
    sol2(2:N,t+1) = Bn\sol2(2:N,t);
end
call_implicit = sol2(S0+1,M+1)
toc



%% change it to be more accurate with smaller ds
X = 60;S0=50;r=0.05;sigma=0.2;Svec=0:0.1:100;tvec=0:0.001:1;oType='c';
tic
deltat = tvec(2)-tvec(1);
ds = Svec(2)-Svec(1);
N = length(Svec)-1;
M = length(tvec)-1;
% j is always follow the S vector pace
j = 0:N;
sig2 = sigma*sigma;
aj = (deltat*j/2).*(r - sig2*j);
bj = 1 + deltat*(sig2*(j.^2) + r);
cj = -(deltat*j/2).*(r + sig2*j);

B = diag(aj(3:N),-1) + diag(bj(2:N)) + diag(cj(2:N-1),1);

[L,U] = lu(B);

sol2=zeros(N+1,M+1);
%set value and set boundary first row to be 0 and last row be maximum value
switch oType
    case 'c'
        sol2(:,1) = max(Svec-X,0);
        for j = 2:M+1
            sol2(N+1,j)=sol2(N+1,j-1)*exp(-r*deltat);
        end
    case 'p'
        sol2(:,1) = max(-Svec+X,0);
        for j = 2:M+1
            sol2(1,j)=sol2(1,j-1)*exp(-r*deltat);
        end
end
        


%detailed loop, get system equalation and solve for next step
offset = zeros(size(B,2),1);
for t = 1:M
    offset(1) = aj(2)*sol2(1,t);
    offset(end) = cj(end)*sol2(end,t);
    sol2(2:N,t+1) = U\(L\(sol2(2:N,t)-offset));
end
call_implicit = sol2(S0/ds+1,M+1)
toc



%% Crank -Nicoloson method

X = 60;S0=50;r=0.05;sigma=0.2;Svec=0:0.1:100;tvec=0:0.0001:1;oType='c';
tic
deltat = tvec(2)-tvec(1);
ds = Svec(2)-Svec(1);
N = length(Svec)-1;
M = length(tvec)-1;
% j is always follow the S vector pace
j = 0:N;
sig2 = sigma*sigma;
dt = deltat
aj = (dt/4)*(sig2*(j.^2) - r*j);
bj = -(dt/2)*(sig2*(j.^2) + r);
cj = (dt/4)*(sig2*(j.^2) + r*j);

D = diag(aj(3:N),-1) + diag(1+bj(2:N)) + diag(cj(2:N-1),1);
C = diag(-aj(3:N),-1) + diag(1-bj(2:N)) + diag(-cj(2:N-1),1);
[L,U] = lu(C);

sol2=zeros(N+1,M+1);
%set value and set boundary first row to be 0 and last row be maximum value
switch oType
    case 'c'
        sol2(:,1) = max(Svec-X,0);
        for j = 2:M+1
            sol2(N+1,j)=sol2(N+1,j-1)*exp(-r*deltat);
        end
    case 'p'
        sol2(:,1) = max(-Svec+X,0);
        for j = 2:M+1
            sol2(1,j)=sol2(1,j-1)*exp(-r*deltat);
        end
end
        


%detailed loop, get system equalation and solve for next step
offset = zeros(size(C,2),1);
for t = 1:M
    offset(1) = aj(2)*sol2(1,t);
    offset(end) = cj(end)*sol2(end,t);
    sol2(2:N,t+1) = U\(L\((D*sol2(2:N,t))+2*offset));
end
call_implicit = sol2(S0/ds+1,M+1)
toc