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

%% Monte Carlo simulation
N_sim = 10000;
tic
[call, se_call] = monte_carlo(S0,K,B,T,r,q,vola_alpha,N_sim);
toc
call
se_call

%% 
N_sim = 20000;
tic
[call, se_call] = monte_carlo(S0,K,B,T,r,q,vola_alpha,N_sim);
toc
call
se_call

%% 
N_sim = 40000;
tic
[call, se_call] = monte_carlo(S0,K,B,T,r,q,vola_alpha,N_sim);
toc
call
se_call

%% 
N_sim = 80000;
tic
[call, se_call] = monte_carlo(S0,K,B,T,r,q,vola_alpha,N_sim);
toc
call
se_call

%% 
N_sim = 160000;
tic
[call, se_call] = monte_carlo(S0,K,B,T,r,q,vola_alpha,N_sim);
toc
call
se_call