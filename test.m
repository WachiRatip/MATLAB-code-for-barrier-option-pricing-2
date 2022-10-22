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

%% FDM: Set the number of grid points
N = 50;         % For the space interval [a,b]
M = 500;        % For the time interval [0,T]

%% solving the Black-Scholes PDE using explicit FDM
[call, V] = explicit(S0,K,B,T,r,q,vola_alpha,N,M);
call

%% solving the Black-Scholes PDE using implicit FDM
[call, V] = implicit(S0,K,B,T,r,q,vola_alpha,N,M);
call

%% solving the Black-Scholes PDE using Crank-Nicolson FDM
[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,M);
call

%% Monte Carlo simulation
N_sim = 10000;  % Number of simulations
[call, se_call] = monte_carlo(S0,K,B,T,r,q,vola_alpha,N_sim);
call