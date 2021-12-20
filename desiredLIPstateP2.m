function [xNorm, U] = desiredLIPstateP2( z0, TS, TD, vdes, targetStepSize, iDomain)
g = 9.81;
lambda = sqrt(g/z0);

% eATs = [ cosh(TS*lambda), sinh(TS*lambda)/lambda;
%          lambda*sinh(TS*lambda), cosh(TS*lambda)];

sigma2 = lambda*tanh(TS/2*lambda);

d = lambda^2*(sech(lambda*TS/2))^2*vdes*(TS + TD)/(lambda^2*TD + 2*sigma2);

X1Fdes = -targetStepSize/(2 + TD*sigma2); %%% left stance
VX1Fdes = sigma2*X1Fdes + d;
[X2Fdes, VX2Fdes, ~] = LIP_sol(lambda, -X1Fdes, VX1Fdes, TS);

if iDomain == 1 %%% right stance%
    xNorm = [X2Fdes; VX2Fdes];
else %%% left stance
    xNorm = [X1Fdes; VX1Fdes];
end

U = 2*xNorm(1) + TD*xNorm(2);
end