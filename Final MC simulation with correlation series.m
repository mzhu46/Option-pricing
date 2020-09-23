r=[0.05;0.05;0.05;0.05;0.05]; sigma=0.2; rho=0.7; T=1; K=100; S0=100;
N = 10^5; % number of MC samples
Omega = eye(5) + rho*(ones(5)-eye(5));
L = chol(Omega,'lower'); % Cholesky factorisation
W = sqrt(T)*L*randn(5,N);
S = S0.*exp((r-0.5*sigma^2)*T + sigma*W);
S1 = 0.2*sum(S,1); % average asset value
F = exp(-r*T)*max(S1-K,0); % call option
val = sum(F)/N; % mean and its std. dev.
sd = sqrt( (sum(F.^2)/N - val.^2)/(N-1) );

corrcoef(S(1,:),S(2,:))


%% with two different mu and sigma

r=[0.03;0.06]; sigma=[0.05;0.1]; rho=0.5; T=1; K=100; S0=[100;120];
N = 10^5; % number of MC samples
Omega = eye(2) + rho*(ones(2)-eye(2));
L = chol(Omega,'lower'); % Cholesky factorisation
W = sqrt(T)*L*randn(2,N);
S = S0.*exp((r-0.5*sigma.*sigma)*T + sigma.*W);
S1 = 0.2*sum(S,1); % average asset value
F = exp(-r*T)*max(S1-K,0); % call option
val = sum(F)/N; % mean and its std. dev.
sd = sqrt( (sum(F.^2)/N - val.^2)/(N-1) );

corrcoef(S(1,:),S(2,:))

%% with multi steps
steps = 50;dt = 1/365;

r=[0.03;0.06]; sigma=[0.05;0.1]; rho=0.5; T=1; K=100; S0=[50;48];

Omega = eye(2) + rho*(ones(2)-eye(2));
L = chol(Omega,'lower'); % Cholesky factorisation
W = sqrt(dt)*L*randn(2,1);

S = zeros ( 2,steps+1);
S(:,1) = S0;

for i = 1 : steps
    W = sqrt(dt)*L*randn(2,1);
    S(:,i+1) = S(:,i).*exp((r-0.5*sigma.*sigma)*dt + sigma.*W);
end

plot (S(1,:))
hold on 
plot (S(2,:))
hold off 
diff= corrcoef(S(1,:),S(2,:))
diff(1,2)
%% verify the coefficient whether equals the rho

Result = zeros(1,1)
for i = 1:10000
    
    steps =100 ;dt = 1/100;
    r=[0.03;0.06]; sigma=[0.05;0.1]; rho=0.5; T=1; K=100; S0=[50;48];

    Omega = eye(2) + rho*(ones(2)-eye(2));
    L = chol(Omega,'lower'); % Cholesky factorisation
    W = sqrt(dt)*L*randn(2,1);

    S = zeros ( 2,steps+1);
    S(:,1) = S0;

    for i = 1 : steps
        W = sqrt(dt)*L*randn(2,1);
        S(:,i+1) = S(:,i).*exp((r-0.5*sigma.*sigma)*dt + sigma.*W);
    end

    diff= corrcoef(S(1,:),S(2,:));
    Result = [Result diff(1,2)];
end
mean(Result)