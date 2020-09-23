%inputs 
X = 60;S0=50;r=0.05;sig=0.2;Svec=0:0.1:100;tvec=0:0.00001:1;
tic
%generate the lattice and boundary condition

M = length(Svec)-1;
N = length(tvec)-1;
dt = tvec(2)-tvec(1);
ds = Svec(2)-Svec(1);
PriceLattice(1:M+1,1:N+1) = nan;
PriceLattice(M+1,:)=0;
%price range is from 0 to 100 and hence 101 nods

for i = 1:M+1
    PriceLattice(i,N+1) = max(0,100/ds+1-i - X/ds);
end

for j = N+1:-1:2
    PriceLattice(1,j-1) = PriceLattice(1,j)*exp(-dt*r);
end
for j = N+1:-1:2
    for i = M: -1:2
        i = 100/ds+1-i;
        a1 = 0.5*dt*(sig^2*i^2-r*i);
        b1 = 1-dt*(sig^2*i^2+r);
        c1=0.5*dt*(sig^2*i^2+r*i);
        i = 100/ds+1-i;
        PriceLattice(i,j-1) = c1*PriceLattice(i-1,j) + b1*PriceLattice(i,j)+...
            a1*PriceLattice(i+1,j);
    end
end

fprintf("Call price with K = %.2f and S0 = %.2f  is %f4" ,X, S0, PriceLattice(S0/ds+1,1)*ds)
toc


%% using matrix form 
%inputs 
X = 60;S0=50;r=0.05;sig=0.2;Svec=0:0.1:100;tvec=0:0.00001:1;
tic

M = length(Svec)-1;
N = length(tvec)-1;
% Get the grid sizes (assuming equi-spaced points)
dt = tvec(2)-tvec(1);
ds = Svec(2)-Svec(1);

% Calculate the coefficients
% To do this we need a vector of j points
j = 1:M-1;
sig2 = sig*sig;
j2 = j.*j;
aj = 0.5*dt*(sig2*j2-r*j);
bj = 1-dt*(sig2*j2+r);
cj = 0.5*dt*(sig2*j2+r*j);

% Pre-allocate the output
price(1:M+1,1:N+1) = nan;

% Specify the boundary conditions

        % Specify the expiry time boundary condition
        price(:,end) = max(Svec-X,0);
        % Put in the minimum and maximum price boundary conditions
        % assuming that the largest value in the Svec is
        % chosen so that the following is true for all time
        price(1,:) = 0;
        price(end,:) = (Svec(end)-X)*exp(-r*tvec(end:-1:1));

% Form the tridiagonal matrix
A = diag(bj);  % Diagonal terms
A(2:M:end) = aj(2:end); % terms below the diagonal
A(M:M:end) = cj(1:end-1); % terms above the diagonal

% Calculate the price at all interior nodes
offsetConstants = [aj(1); cj(end)];
for i = N:-1:1
    price(2:end-1,i) = A*price(2:end-1,i+1);
    % Offset the first and last terms
    price([2 end-1],i) = price([2 end-1],i) + ...
        offsetConstants.*price([1 end],i+1);
end

% Calculate the option price
fprintf("Call price with K = %.2f and S0 = %.2f  is %f4" ,X, S0, price(S0/ds+1,1))
toc

