function [stepL, dstepL] = P1stab(z0, xNow, dxNow, Xdes, Udes, TS, TD) 
%%
%%%% xNow: current COM pos (sagittal) relative to stance foot, 
%%%% dxNow: velocity of the COM
%%%% TS, TD: duration of the SSP and DSP
%%%% stepL: desired target step length
%%%% vdes: desired velocity
%%%% z0 appro height of the COM 
g = 9.81; 
lambda = sqrt(g/z0); 
X = [xNow; dxNow]; 
% %%%%%% different gains %%%%%%%%%%%%%%%
% % %% K deadbeat 
k1 = 1;
k2 = TD + coth(TS*lambda)/lambda;
K = [k1 , ...
     k2 ];

% %% LQR
% N = zeros(2,1);
% R = 1;
% Q1 = 1e+3; 
% Q2 = 1e+2; 
% Q = [Q1, 0; 0, Q2];
% [K, ~, ~] = dlqr(Ahat, Bhat, Q, R, N); %%% 
% K = -K; %%%% sign convention is different here. 
%%
stepL = K*(X - Xdes) + Udes;
dX = [dxNow;g/z0*xNow ];
dstepL = K*dX; 
end 