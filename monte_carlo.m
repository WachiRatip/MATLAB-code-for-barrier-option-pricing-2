% A function to calculate the price of a barrier call option price
% for a local volatility Black-Scholes model
% using the Monte Carlo simulation
%
% monte_carlo(spot,strike,barrier,term,r,d,vola_alpha,N)
%
% Inputs: stock     - stock price
%       : strike    - strike price
%       : barrier   - barrier level
%       : term      - time to maturity (years)
%       : r         - interest rate (%)	
%       : d         - dividend yield (%)
%       : vola_alpha- the local volatility surface parameter
%       : N         - number of simulations
%
% Outputs: call     - The call option price corresponding to the stock price
%        : se_call  - The standard error of the option prices 
%                     corresponding to the stock price


%% implements the Monte Carlo simulation
function [call, se_call] = monte_carlo(stock,strike,barrier,term, ...
    r,d,vola_alpha,N)
% Rearrange parameters; change from % to decimal
r = r/100;
d = d/100;

% compute sigma
sigma = 0.25*exp( -term ) * ( (100/stock)^(vola_alpha) );

% placeholder to store sample of stock prices at maturity
S = zeros(N,2);
% placeholder to store sample of discounted call payoff
X = zeros(N,2);

% Monte-Carlo scheme
for k=1:N
    % start sampling Z ~ N(0,1)
    Z = randn;
    
    % calculate the stock price at time T
    S(k,1) = stock*exp( (r-d-0.5*sigma^2)*term + sigma*sqrt(term)*Z );
    % compute the call option price;
    % if S(k) > B, then X(k) = 0 otherwise exp(-r*T)*max((S(k)-K),0)
    X(k,1) = exp(-r*term) * max((S(k,1)-strike), 0) * (S(k,1)<=barrier);

    % antithetic variates
    S(k,2) = stock*exp( (r-d-0.5*sigma^2)*term - sigma*sqrt(term)*Z );
    X(k,2) = exp(-r*term) * max((S(k,2)-strike), 0) * (S(k,2)<=barrier);

end
% compute the average call & the standard error
call = mean(mean(X, 2));
se_call = sqrt( (sum(mean(X.^2, 2))-(N*call^2)) / (N*(N-1)) );

end