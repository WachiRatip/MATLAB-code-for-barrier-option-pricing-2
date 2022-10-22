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

%% solving the Black-Scholes PDE using explicit FDM
N = 50;
[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,100);
call

[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,200);
call

[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,300);
call

[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,400);
call

[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,500);
call

%% 
N = 100;
[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,100);
call

[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,200);
call

[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,300);
call

[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,400);
call

[call, V] = crank(S0,K,B,T,r,q,vola_alpha,N,500);
call