function [xNorm, U] = desiredLIPstateP1( z0, TS, TD, vdes)
g = 9.81; 
lambda = sqrt(g/z0); 
            
eATs = [ cosh(TS*lambda), sinh(TS*lambda)/lambda;
    lambda*sinh(TS*lambda), cosh(TS*lambda)];

sigma1 = lambda*coth(TS/2*lambda);
Bhat = eATs*[-1;
    0];
Ahat = eATs*[1, TD;
            0, 1];
pF = vdes/sigma1;
vF = vdes; 
xNorm = [pF; vF];
U = -transpose((Ahat - eye(2))*xNorm)*[1;0]/Bhat(1);

end