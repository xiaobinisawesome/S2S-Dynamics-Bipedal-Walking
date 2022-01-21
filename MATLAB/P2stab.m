function [stepW, dstepW] = P2stab(z0, yNow, dyNow, Ydes, Udes, TS, TD) 

%%%% yNow: current COM pos (coronal) relative to stance foot, 
%%%% dyNow: velocity of the COM
%%%% TS, TD: duration of the SSP and DSP
%%%% stepW: desired target step width
%%%% vdes: desired velocity
%%%% z0 appro height of the COM 

g = 9.81; 
lambda = sqrt(g/z0); 
k1 = 1;
k2 = TD + coth(TS*lambda)/lambda; %%
K = [k1 , ...
      k2 ];
ynow = [yNow; dyNow]; 
stepW = K*(ynow - Ydes) + Udes;

dY = [dyNow; g/z0*yNow];
dstepW = K*dY; 
end