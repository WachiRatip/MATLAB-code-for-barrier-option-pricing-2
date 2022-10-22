%% Clear memory and colsone output
clc
clear

%% the problem parameters 
S0 = 100;           % spot price (in British Pound)
K = 90;             % strike price (in British Pound)
B = 130;            % barrier level (in British Pound)
r = 3;              % risk-free rate (in %)
q = 5;              % dividend yield (in %)
T = 0.5;            % time to maturity (years)
vola_alpha = 0.35;  % the local volatility alpha

%% the model parameters
% FDM: Set the number of grid points
N = 50;         % For the space interval [a,b]
M = 500;        % For the time interval [0,T]

% Monte Carlo: Set the number of simulations
N_sim = 10000;  % Number of simulations
M_sample = 100; % Number of discrete time steps

%% solving the Black-Scholes PDE using explicit FDM
[call, V] = explicit(S0,K,B,T,r,q,vola_alpha,N,M);

% plot of FDM solution
subplot(2,1,1)
dS = B/N;
S = (0+dS:dS:B-dS)';
plot(S,V(:,M+1),'LineWidth',2)
hold on
plot(S0,call,'r*') % mark a call option price on the plot 
title('European Call price, Explicit - BS formula')
xlabel('Stock price')
ylabel('Call price')

subplot(2,1,2)
t = T-(0:T/M:T);
[X,Y] = meshgrid(t,S);
surf(X,Y,V,'LineStyle','none');
xlabel('t')
ylabel('S(t)')
zlabel('V(S,t)')

%% solving the Black-Scholes PDE using implicit FDM
[call, V] = implicit(S0,K,B,T,r,q,vola_alpha,N,M);

% plot of FDM solution
subplot(2,1,1)
dS = B/N;
S = (0+dS:dS:B-dS)';
plot(S,V(:,M+1),'LineWidth',2)
hold on
plot(S0,call,'r*') % mark a call option price on the plot 
title('European Call price, Implicit - BS formula')
xlabel('Stock price')
ylabel('Call price')

subplot(2,1,2)
t = T-(0:T/M:T);
[X,Y] = meshgrid(t,S);
surf(X,Y,V,'LineStyle','none');
xlabel('t')
ylabel('S(t)')
zlabel('V(S,t)')

%% solving the Black-Scholes PDE using Crank-Nicolson FDM
[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,M);

% plot of FDM solution
subplot(2,1,1)
dS = B/N;
S = (0+dS:dS:B-dS)';
plot(S,V(:,M+1),'LineWidth',2)
hold on
plot(S0,call,'r*') % mark a call option price on the plot 
title('European Call price, Crank-Nicolson - BS formula')
xlabel('Stock price')
ylabel('Call price')

subplot(2,1,2)
t = T-(0:T/M:T);
[X,Y] = meshgrid(t,S);
surf(X,Y,V,'LineStyle','none');
xlabel('t')
ylabel('S(t)')
zlabel('V(S,t)')

%% Monte Carlo simulation
[call, se_call] = monte_carlo2(S0,K,B,T,r,q,vola_alpha,N_sim, M_sample);

% compute the 95p confidance interval
ci_low = call - 1.96*se_call;
ci_high = call + 1.96*se_call;