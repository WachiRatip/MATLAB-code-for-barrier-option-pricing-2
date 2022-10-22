% A function to calculate the price of a barrier call option price
% for a local volatility Black-Scholes model
% using the Monte Carlo simulation
%
% monte_carlo(spot,strike,barrier,term,r,d,vola_alpha,N,M)
%
% Inputs: stock     - stock price
%       : strike    - strike price
%       : barrier   - barrier level
%       : term      - time to maturity (years)
%       : r         - interest rate (%)	
%       : d         - dividend yield (%)
%       : vola_alpha- the local volatility surface parameter
%       : N         - number of simulations
%       : M         - number of discrete time steps
%
% Outputs: call     - The call option price corresponding to the stock price
%        : se_call  - The standard error of the option prices 
%                     corresponding to the stock price


%% implements the Monte Carlo simulation
function [call, se_call] = monte_carlo2(stock,strike,barrier,term, ...
    r,d,vola_alpha,N,M)
% Rearrange parameters; change from % to decimal
r = r/100;
d = d/100;

% Time interval (0 <= tau <= T)
dtau = term/M;
tau = 0:dtau:term;

% placeholder to store sample of stock prices at maturity
S = zeros(M+1,2);
% placeholder to store sample of discounted call payoff
X = zeros(N,2);

% Monte-Carlo scheme
for i=1:N
    S(1,1) = stock;
    S(1,2) = stock;
    for k=1:M
        % start sampling Z ~ N(0,1)
        Z = randn;
        % calculate the local volatility surface function; sigma(S, t).
        sigma = 0.25*exp( -tau(k) ) * ( (100/S(k,1))^(vola_alpha) );
        % calculate the stock price at time T
        S(k+1,1) = S(k,1)*exp( (r-d-0.5*sigma^2)*dtau + sigma*sqrt(dtau)*Z );
        % antithetic variates
        sigma = 0.25*exp( -tau(k) ) * ( (100/S(k,2))^(vola_alpha) );
        S(k+1,2) = S(k,2)*exp( (r-d-0.5*sigma^2)*dtau - sigma*sqrt(dtau)*Z );
    end
    
    % compute the call option price;
    % if S(T) > B, then X(i) = 0 otherwise exp(-r*T)*max((S(T)-K),0)
    X(i,1) = exp(-r*term) * max((S(M+1,1)-strike), 0) * (S(M+1,1)<=barrier);
    X(i,2) = exp(-r*term) * max((S(M+1,2)-strike), 0) * (S(M+1,2)<=barrier);

end
% compute the average call & the standard error
call = mean(mean(X, 2));
se_call = sqrt( (sum(mean(X.^2, 2))-(N*call^2)) / (N*(N-1)) );
plot(S(:,:))
end