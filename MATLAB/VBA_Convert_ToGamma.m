function [a, b] = VBA_Convert_ToGamma(mu, sig)
%[ a,b ] = VBA_Convert_ToGamma( mu,sig )
%
% Converts mu and sig defining a prior distribution over noise standardard
% deviation s, to a gamma distriribution over the associated precision p
% defined by a and b with p = 1/s^2

% IN:
%   mu,sig: mean and standard deviation of prior distribution over s
% OUT:
%   a,b: Can either be a_sigma/b_sigma or a_alpha/b_alpha defining
%   prior distributions over measurement/system noise precision
% 
%  M. Eichenlaub 04/12/2019

if sig == 0     % Return Inf and 0 if stochstic inversion shall be disabled
    a = Inf;
    b = 0;
else

% Check mu and sig to prevent errors in numerical approximation        
assert(nargin==2,' *** VBA_Convert_ToGamma: Both mu and sigma must be provided')
assert(mu<2e-3 || mu>1e4,' *** VBA_Convert_ToGamma: mu must be greater than 2e-3 and smaller than 1e4')
assert(sig/mu<2e-3 || sig/mu>50,' *** VBA_Convert_ToGamma: sig must be greater than 2e-3*mu and smaller than 50*mu')

% Upper bound for a
a0 = 1/8*(1+sqrt(49+mu^4/sig^4+50*mu^2/sig^2)+mu^2/sig^2);

% Find a by minimizing the function D for 1<a<a0
opt = optimset('MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-5);
a = fminbnd(@(x)D(x,mu,sig),1,a0,opt);

% Calculate b
b = (mu./(exp(gammaln(a-0.5)-gammaln(a)))).^2;

end

end

function [ D ] = D(a,mu,sig)
% Function has to be minimized with respect to a

S = exp(2*(gammaln(a-0.5)-gammaln(a)));
D = mu^2./S - sig^2./((1./(a-1))-S);
D = log(D.^2+1);

end