% A function to calculate the price of a barrier call option price
% for a local volatility Black-Scholes model
% using the Crank-Nicolson finite difference method
%
% crank(spot,strike,barrier,term,r,d,vola_alpha,N,M)
%
% Inputs: stock     - stock price
%       : strike    - strike price
%       : barrier   - barrier level
%       : term      - time to maturity (years)
%       : r         - interest rate (%)	
%       : d         - dividend yield (%)
%       : vola_alpha- the local volatility surface parameter
%       : N         - number of grid points for the space interval [a,b]
%       : M         - number of grid points for the time interval [0,T]
%
% Outputs: call     - The call option price corresponding to the stock price
%        : V        - The option prices calculated from BS via FDM.


%% implements the Crank-Nicolson scheme

function [call, V] = crank(stock,strike,barrier,term,r,d,vola_alpha,N,M)
% Rearrange parameters; change from % to decimal
r = r/100;
d = d/100;

% Grid points construction
% Space interval (Smin <= S <= Smax); the minimal and maximal stock prices
Smin = 0;           % a in [a,b]
Smax = barrier;     % b in [a,b]
dS = (Smax-Smin)/N; % (b-a)/N
S = (Smin+dS:dS:Smax-dS)';

% Time interval (0 <= tau <= T)
dtau = term/M;
tau = 0:dtau:term;

% Calculate the local volatility surface function; sigma(S, t).
sigma = zeros(N-1,M+1);
for k=1:M+1
    % the local volatility surface: 
    % sigma(S,t) = 0.25*e^{-t}*(100/S)^{alpha}
    sigma(:,k) = 0.25*exp( -tau(:,k) ) * ( (100./S).^(vola_alpha) );
end

% Calculate the alpha and beta paramters for V(S, tau)
alpha = zeros(N-1,M+1);
beta = zeros(N-1,M+1);
for k=1:M+1
    % alpha = (sigma^2*S^2/2)*(dt/dS^2)
    alpha(:,k) = (0.5 * (sigma(:,k).^2) .* (S.^2) * dtau) / (dS^2);
    % beta = (r-d)S/2*(dt/dS)
    beta(:,k) = ((r - d) * S * dtau) / (2*dS);
end

% the explicit finite difference scheme's
% lower, main and upper diagonal components of the tridiagonal matrix 
le = alpha - beta;
de = 2 - r*dtau - 2*alpha;
ue = alpha + beta;

% the implicit finite difference scheme's
% lower, main and upper diagonal components of the tridiagonal matrix 
li = beta - alpha;
di = 2 + r*dtau + 2*alpha;
ui = - alpha - beta;

% Solve BS PDE by using the Crank-Nicolson scheme
V = zeros(N-1,M+1);
% a solution matrix with initial condition
V(:,1) = max(S-strike, 0);

for k=1:M
    % Construct AE: the explicit tridiagonal matrix, (AE+I).
    AE = diag(de(:,k)) + diag(ue(1:N-2,k),1) + diag(le(2:N-1,k),-1);

    % Construct AI: the implicit tridiagonal matrix, (AI+I).
    AI = diag(di(:,k+1)) + diag(ui(1:N-2,k+1),1) + diag(li(2:N-1,k+1),-1);

    % Iteration to compute the value V(j,k)
    % with zero boundary conditions, i.e., V(Smin) = V(Smax) = 0
    V(:,k+1) = AI\(AE*V(:,k));
end % k

% Interpolation to find the call price at the stock price: S0 
call = interp1(S,V(:,M+1),stock);

end