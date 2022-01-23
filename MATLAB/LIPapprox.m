classdef LIPapprox < handle
    %%%% ---- Xiaobin Xiong ----------
    %%%% define a LIP approximation class
    %%%% which is used for aSLIP walking, Cassie walking
    %%%% many ideas in this file. 
    %%%% updated 1/23/2022, I may have time to do a clean up later. 
    
    properties
        g = 9.81;
        mass = 33.3; 
        TS % duration of SSP
        TD % duration of DSP
        z % nominal height of COM
        eAT
        theta = 0 %%% local slope
        vxdes = 0 %%% desired net velocity
        lambda
        
        TSavg  %%% averaged duration during walking 
        TDavg 
        
        NominalTraj = struct('t', [], 'x', [], 'x2f', [], 'dx', []); %%% 
        
        %%%% control type 
        goalBehavior%%% = 'localBehavior'  %%/ 'globalBehavior' 
        whichController %%% 'Deadbeat' 'MPC'
        whichGain = 'deadbeat' %% default
        stanceLeg
        stanceLegNum
        
        gaitType %%% P1
        LIPacc =0
        targetStepSize = 0; %%% of P2 orbit
    end
    
    properties % periodic behaviors
        
        %%% orbital slopes
        sigma1   %% sigma1   
        sigma2   %% sigma2
                
        Xdes = zeros(2,1) % boundary state for P1 orbits
        Udes = 0 ; 
        
        x1Fdes = 0 % selected boundary position for P2 orbits. 
        pLeftDes = -0.05; %%% same as x1Fdes, todo: clear this up
        XLdes = zeros(2,1);  %%% for P2 orbits
        XRdes = zeros(2,1);
        d2 = 0; 
        LIP_dxf_previous = 0;
        
        % struct for logging 
        HLIP = struct('t', [], 'x', [], 'x2f',[], 'dx', [], 'pred_p', [],'pred_v', [], 'u', []);
        robotXF = struct('t',[], 'p', [], 'x', [], 'v',[], 'u', [], 'w', []);
        
        X0previous = [];
        
        LIPkd = 0;   %%%%%% feedback raibert_kd  %%% add derivative term for comparison
        Tsol
        Xsol
        Ahat0 %%% original LIP, constant
        Bhat0
        
        %%% used for global position control
        A3
        B3
        
        %%% used for learning
        Anow
        Bnow
        
        %%% linear regressed S2S
        Abar = zeros(2,2);
        Bbar = zeros(2,1); 
        Cbar = zeros(2,1); 
        episonMax = zeros(2,1)
        pushDisturbance = zeros(2,1); 
        ifRS2S = 0; 
        
        
        cvd %%%% the offset of input (stepsize) to work on the error dynamics 
        %%%  robust MPC setup
        Xset %%% state set
        Uset %%% input set
        Wset %%% disturbance set, model different 
        Eset %%% disturbance invariant set 
        X3set %%% include global position. 
        W3set 
        E3set 
        
        %%%%% discrete LIP states %%%%%%%%%%
        Xint = 0;
        XLIP = zeros(2,1);
        XLIPlist = [];  % saved
        X3LIP = zeros(3,1);
        X3LIPlist = [];
    end
    
    properties % controllers
        LIP2SLIP_deltaVec = [];   
        %%%% 
        K
        K3
        %%%% deadbeat
        Kdeadbeat 
        K3deadbeat
        %%%% LQR
        KLQR = zeros(1,2); 
        K3LQR = zeros(1,3); 
        %%% adaptive
        regreeEst = zeros(2);
        C = zeros(2,1);
        
        %%%%%%%%%% system level synthesis %%%%%%%%%
        N_fir  = 5;%%%%%%%%%% finite impluse response
        w_est_vec = []
        SLS_XrefGain = [];
        SLS_Ugain = [];
        SLS_ulog =[];
        x_ref_k = zeros(2,1); 
        x_ref_km1 = zeros(2,1); 
        %%% robust
        % MPC based controllers
        MPC_Lresults = [] % save for the MPC step length vector
        targetStepLength = 0 % target step length for the current step
        targetStepLengthVec
        targetStepLengthOptVec %%% optimized vecl
        MPCsol = struct('xsol', [], 'usol', [])        
 
        %%% optimization based 
        optStepping 
        
        %%% 
        latestStepL %%%%  most recent calcuation. 
    end
    
    properties %%%% parameter ID, gradient Method %%%%%%
        Zid
        EpsilonID
        XID
        Am 
        Gamma 
        ThetaID 
        ThetaIDlist 
        X_tilde = zeros(2, 1); 
        Klist  %%% what is the K here.  
        Know  %% current state-feedback gain
        
        Emeasurement = [0.0;0]; %zeros(2,1); %%% measurement error on the LIP state. %%%%% used for testing robustness
    end 
    
    methods
        
        function obj = LIPapprox(TS, TD, z)
            %%%%%
            obj.TS = TS;
            obj.TD = TD;
            obj.z = z;
            
            obj.lambda = sqrt(obj.g/obj.z); 
            
            eATs = [ cosh(TS*obj.lambda), sinh(TS*obj.lambda)/obj.lambda;
                obj.lambda*sinh(TS*obj.lambda), cosh(TS*obj.lambda)];
            
            %%%% step-to-step dynamics %%%%%%
            B_hat = eATs*[-1;
                0];
            A_hat = eATs*[1, TD;
                0, 1];  
            
            obj.Ahat0 = A_hat; 
            obj.Bhat0 = B_hat;
            
            %%%%
         %   obj.optStepping = optimizationStepping(0);
            
        %    obj.initializeLIP(TS, TD)
        %    obj.initializeParameterID();
        %    obj.initializeSet();
            %%% 
            obj.w_est_vec = zeros(2,obj.N_fir); 
        end
        
        function obj = initializeLIP(obj) 
            eATs = [ cosh( obj.TS*obj.lambda), sinh( obj.TS*obj.lambda)/obj.lambda;
                obj.lambda*sinh( obj.TS*obj.lambda), cosh( obj.TS*obj.lambda)];
            %%% 
            %   obj.optStepping.updateLIP(obj.lambda, obj.TS, obj.TD)
            obj.Bhat0 = eATs*[-1;
                              0];
            obj.Ahat0 = eATs*[1,  obj.TD;
                              0, 1];
            
            obj.A3 = [1, obj.Ahat0(1,1)-1, obj.Ahat0(1,2);
                zeros(2,1), obj.Ahat0];
            obj.B3 = [obj.Bhat0(1)+1;
                obj.Bhat0];
            
            obj.sigma2 = obj.lambda*tanh(obj.TS/2*obj.lambda);
            obj.sigma1 = obj.lambda*coth(obj.TS/2*obj.lambda);
            
            %%%%% for P2 orbits %%%%%%%%%
            % characteristic line dotX = k x + d;
            d = obj.lambda^2*(sech(obj.lambda*obj.TS/2))^2*obj.vxdes*(obj.TS+obj.TD)/(obj.lambda^2*obj.TD + 2*obj.sigma2); 
            obj.d2 = d; 
            % select the boundary condition. % x1_F 
            pLdes = obj.pLeftDes ; %-0.05; % m  custom selection 
            vLdes = obj.sigma2*pLdes + d; 
            [pRdes, vRdes, ~] = LIP_sol(obj.lambda, -pLdes, vLdes, obj.TS); 
            obj.XLdes = [pLdes; vLdes]; 
            obj.XRdes = [pRdes; vRdes]; 
           
            %%% update dLQR gain
            %%%% LQR %%%%%%%%%%%%%%%%%%
            N = zeros(2,1);
            R = 4;
            Q = blkdiag(1,0.5);
             %Q = blkdiag(0.01,0.01);
             
            [obj.KLQR, ~, ~] = dlqr(obj.Ahat0, obj.Bhat0, Q, R, N);
            obj.KLQR = -obj.KLQR; %%% sign convention
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N = zeros(3,1);
            R = 1; %1e+2;
           % Q = 1e+1*blkdiag(1,.1,1);
            Q = 1e+1*blkdiag(10,.1,1);
            [obj.K3LQR, ~, ~] = dlqr(obj.A3, obj.B3, Q, R, N);
            obj.K3LQR = - obj.K3LQR;
            %%%%%%%%%%%%%%%%%
            obj.Kdeadbeat = [1, obj.TD + coth(obj.TS*obj.lambda)/obj.lambda] ;
             
            k1 = 1/(-2 + 2*cosh(obj.TS*obj.lambda) + obj.TD*obj.lambda*sinh(obj.TS*obj.lambda));
            k2 = 1; 
            k3 = (4*obj.TD*obj.lambda + coth(obj.TS*obj.lambda/2)*(3+2*obj.TD^2*obj.lambda^2+2*obj.TD*obj.lambda*coth(obj.TS*obj.lambda/2)) + tanh(obj.TS*obj.lambda/2))/...
                (2*obj.lambda*(2 + obj.TD*obj.lambda*coth(obj.TS*obj.lambda/2)));
            obj.K3deadbeat = [k1, k2, k3];  %%% robustifying
            
            switch obj.whichGain 
                case 'deadbeat'
                    obj.K = obj.Kdeadbeat;                     
                    obj.K3 = obj.K3deadbeat; 
                case 'LQR'
                    obj.K = obj.KLQR;
                  	obj.K3 = obj.K3LQR; 
            end 
        end 
        
        function initializeSet(obj)
            w1 = 0.05;  w2 = 0.1; %%% model difference
            xlim = [0.5, 1];
            ulim = 0.6;
            %%%% set definition. %%%%
            W = Polyhedron([w1, w2; -w1, w2; w1, -w2; -w1, -w2]); W.computeHRep;
            X = Polyhedron([xlim(1), xlim(2); -xlim(1), xlim(2); xlim(1), -xlim(2); -xlim(1), -xlim(2)]); X.computeHRep;
            U = Polyhedron([ulim; -ulim]); U.computeHRep;
            %%%
 
            Acl = obj.Ahat0 + obj.Bhat0*obj.Kdeadbeat;
            E = computeInvariant(Acl, W);
            %%%
            nn = 100; %%% some big number %%%
            X3 = Polyhedron([nn, xlim(1), xlim(2); ...
                nn, -xlim(1), xlim(2); ...
                nn, xlim(1), -xlim(2); ...
                nn, -xlim(1), -xlim(2);...
                -nn, xlim(1), xlim(2); ...
                -nn, -xlim(1), xlim(2); ...
                -nn, xlim(1), -xlim(2); ...
                -nn, -xlim(1), -xlim(2)]);
            X3.computeHRep;
            w0vec = [ obj.A3(1,2)*w1 + obj.A3(1,3)*w2;
                obj.A3(1,2)*w1 - obj.A3(1,3)*w2;
                -obj.A3(1,2)*w1 + obj.A3(1,3)*w2;
                -obj.A3(1,2)*w1 - obj.A3(1,3)*w2;
                ];  %%%
            w0 = max(abs(w0vec));
            W3 = Polyhedron([w0, w1, w2; ...
                             w0, -w1, w2; ...
                             w0, w1, -w2; ...
                             w0, -w1, -w2;...
                             -w0, w1, w2; ...
                             -w0, -w1, w2; ...
                             -w0, w1, -w2; ...
                             -w0, -w1, -w2]); W3.computeHRep;
            Acl3 = obj.A3 + obj.B3*obj.K3deadbeat; 
          %  Acl3 = obj.A3 + obj.B3*obj.K3LQR; 
            E3 = computeInvariant(Acl3, W3); 
           
            obj.Xset = X; 
            obj.Eset = E; 
            obj.Wset = W; 
            
            obj.X3set = X3; 
            obj.E3set = E3; 
            obj.W3set = W3; 
            
            obj.Uset = U; 
        end 
        
        function setDesiredVelocity(obj, v) 
            %%% set desired velocity
            obj.vxdes = v;
            % above is averaged velocity. The SSP xF_dot is the following.
            DistSum = obj.vxdes*(obj.TS + obj.TD); % total traveled distance in SSP and DSP.
            
            switch obj.gaitType
                case 'P1'
                    obj.vxdes = DistSum/(2/obj.sigma1 + obj.TD);
                    obj.Xdes = [obj.vxdes/obj.sigma1; obj.vxdes];
                    obj.Udes = obj.Xdes(1)*2 + obj.Xdes(2)*obj.TD; 
                case 'P2'
                    obj.d2 = obj.lambda^2*(sech(obj.lambda*obj.TS/2))^2*obj.vxdes*(obj.TS+obj.TD)...
                        /(obj.lambda^2*obj.TD + 2*obj.sigma2);
                    X1Fdes = obj.targetStepSize/(2 + obj.TD*obj.sigma2);
                    obj.x1Fdes = X1Fdes;
                    VX1Fdes = obj.sigma2*X1Fdes + obj.d2;
                    obj.XLdes = [obj.x1Fdes; VX1Fdes];
                    [x2Fdes, VX2Fdes, ~] = LIP_sol(obj.lambda, -X1Fdes, VX1Fdes, obj.TS);
                    obj.XRdes = [x2Fdes; VX2Fdes];
            end
            
            if obj.ifRS2S == 1 
                setDesiredVelocityRS2S(v)
            end
        end 
        
        function setDesiredVelocityRS2S(obj, v)
            %%% set desired velocity
            obj.vxdes = v;
            % above is averaged velocity. The SSP xF_dot is the following.
            DistSum = obj.vxdes*(obj.TS + obj.TD); % total traveled distance in SSP and DSP.
            
            switch obj.gaitType
                case 'P1'
                     obj.Udes = DistSum;
                    obj.Xdes = inv(eye(2) - obj.Abar)*(obj.Bbar*obj.udes + obj.Cbar); 
                case 'P2'
                    obj.uLdes = obj.targetStepSize; 
                    obj.uRdes = 2*DistSum - uL; 
                    obj.XLdes = inv(eye(2) - obj.Abar*obj.Abar)*((obj.Abar*obj.Bbar - obj.Bbar)*obj.uLdes + obj.Bbar*2*vy*obj.T + (obj.Abar + eye(2))*obj.Cbar); 
                    obj.XRdes = inv(eye(2) - obj.Abar*obj.Abar)*((obj.Abar*obj.Bbar - obj.Bbar)*obj.uRdes + obj.Bbar*2*vy*obj.T + (obj.Abar + eye(2))*obj.Cbar); 
             end
        end
 
        function prepareSLS(obj, SLS)
            obj.Abar = SLS.Abar;
            obj.Bbar = SLS.Bbar;
            obj.Cbar = SLS.Cbar;
            obj.episonMax = SLS.Wm;
            obj.Xdes = SLS.Xdes;
            obj.Udes = SLS.udes;
            
            a = obj.Abar(1,1); b = obj.Abar(1,2); c = obj.Abar(2,1); d = obj.Abar(2,2);
            e = obj.Bbar(1); f = obj.Bbar(2);
            k1 = (a*c*e + c*d*e - a^2*f - b*c*f)/(-c*e^2 + a*e*f - d*e*f + b*f^2);
            k2 = (b*c*e + d^2*e - a*b*f - b*d*f)/(-c*e^2 + a*e*f - d*e*f + b*f^2);
            
            obj.Kdeadbeat = [k1, k2];
        end 
        
        function setZ(obj, z) 
            obj.z = z; 
        end 
        
        function setP2orbitUleft(obj, u)
            obj.targetStepSize = u; 
            obj.pLeftDes = (u - obj.TD*obj.d2)/(2+ obj.TD*obj.sigma2); 
            obj.initializeLIP( ); 
        end
        
        function setP2orbitPLeft(obj, p) 
            obj.pLeftDes = p; 
            obj.initializeLIP( ); 
        end 
        
        function getPushDisturbance(obj, F)
            %%%% get push disturbance on H_LIP, push for SSP of H-LIP
            obj.pushDisturbance = F/obj.mass/obj.lambda*sinh(obj.TS*obj.lambda)*[1/obj.sigma1;1];
        end
        
        function [x, x2f, dx] = LIPbasedEst(obj, x, x2f,dx, t) 
            %%% est the preimpact state
            %%%% x2f: com pos related to stance pivot 
            %%% x: global pos 
            %%% dx: velocity
             eAt = [cosh(t*obj.lambda), sinh(t*obj.lambda)/obj.lambda;
                obj.lambda*sinh(t*obj.lambda), cosh(t*obj.lambda)];
             stanceP = x - x2f; 
             XF = eAt*[x2f;dx]; 
             x2f = XF(1); 
             dx = XF(2); 
             x = stanceP + x2f;
        end 
        
    end
    
    methods %%%%% if adaptive control, parameter ID
        
        function initializeParameterID(obj)
            %%%%%%%%%%%% Parameter ID on the linear system from actual walking data %%%%%%%%%%%%
            obj.Am = blkdiag(0.05, 0.05);
            obj.Gamma = blkdiag(2, 2); %  gain
            obj.ThetaID = transpose([obj.Ahat0 - obj.Am, obj.Bhat0]);
            obj.ThetaIDlist = reshape(obj.ThetaID', 1, []);  % initial Theta0
            obj.Klist = [];
            obj.Know = [];
        end
        
        function updateParameterID(obj)
            %% update parameters at SSP-
            if length(obj.LIP.xf_list) > 1
                Xpre = [obj.LIP.xf_list(end-1); obj.LIP.dxf_list(end-1)] + obj.Emeasurement;
                previousL = obj.stepLengthSequence(end-1);
            else
                Xpre = [0;0];
                previousL = 0;
            end
            Xnow = [obj.LIP.xf_list(end); obj.LIP.dxf_list(end)] + obj.Emeasurement;
            
            phiNow = [Xpre; previousL];  %%%%% recheck the sequence here!!!
            ThetaPre = obj.ThetaID;
            obj.Zid = Xnow - obj.Am *Xpre;
            obj.EpsilonID = obj.Zid - ThetaPre'*phiNow;
            obj.ThetaID = ThetaPre + transpose(obj.Gamma*obj.EpsilonID*phiNow');
            obj.ThetaIDlist = [obj.ThetaIDlist; reshape(obj.ThetaID' , 1, [])];
            obj.Klist = [obj.Klist, obj.Know];
            
            %%%%%
            %             phiNow = [Xnow; obj.curstepLength];  %%%%% recheck the sequence here!!!
            %             ThetaPre = obj.ThetaID;
            %             obj.Zid = Xnow - obj.Am *Xpre;
            %             obj.EpsilonID = obj.Zid - ThetaPre'*phiNow;
            %             obj.ThetaID = ThetaPre + transpose(obj.Gamma*obj.EpsilonID*phiNow');
            %             obj.ThetaIDlist = [obj.ThetaIDlist; reshape(obj.ThetaID' , 1, [])];
            %             obj.Klist = [obj.Klist, obj.Know];
            %
            %%%% temp try
            %             thetaNow = transpose(obj.ThetaID);
            %             Anow = thetaNow(:, 1:2) + obj.Am;
            %             Bnow = thetaNow(:, 3);
            %             if length(obj.LIP.xf_list) > 1
            %                 Xpre = [obj.LIP.xf_list(end-1); obj.LIP.dxf_list(end-1)] + obj.Emeasurement;
            %             else
            %                 Xpre = [0;0] + obj.Emeasurement;
            %             end
            %             Xnow = [obj.LIP.xf_list(end); obj.LIP.dxf_list(end)]+ obj.Emeasurement;
            %             obj.delta = inv(eye(2) - Anow)*(Xnow - Anow*Xpre - Bnow*obj.curstepLength);
            %%%%%% with parameter update %%%%%
            thetaNow = transpose(obj.ThetaID);
            obj.Anow = thetaNow(:, 1:2) + obj.Am;
            obj.Bnow = thetaNow(:, 3);
        end
        
        function updateMeasurementID(obj)
            %% update parameters at SSP-
            if length(obj.LIP.xf_list) > 1
                Xpre = [obj.LIP.xf_list(end-1); obj.LIP.dxf_list(end-1)] + obj.Emeasurement;
                previousL = obj.stepLengthSequence(end-1);
                if length(obj.LIP.xf_list) > 2
                    prepreL = obj.stepLengthSequence(end-2);
                    Xprepre = [obj.LIP.xf_list(end-2); obj.LIP.dxf_list(end-2)] + obj.Emeasurement;
                else
                    prepreL = 0;
                    Xprepre = zeros(2,1);
                end
            else
                Xpre = [0;0] + obj.Emeasurement;
                Xprepre = [0;0] + obj.Emeasurement;
                previousL = 0;
                prepreL = 0;
            end
            Xnow = [obj.LIP.xf_list(end); obj.LIP.dxf_list(end)]+ obj.Emeasurement;
            dynamicsDelta = Xnow - obj.Ahat0*Xpre - obj.Bhat0*previousL;
            obj.deltaVec = [obj.deltaVec, dynamicsDelta];
            %  obj.delta = 0.5*(obj.delta + inv(eye(2) - obj.Ahat0)*(Xnow - obj.Ahat0*Xpre - obj.Bhat0*previousL));
            
            %%%%
            Kest = -inv(eye(2) - obj.Ahat0)*obj.Ahat0;
            obj.X_tilde = obj.Ahat0*Xpre + obj.Bhat0*previousL + Kest*(Xpre - obj.Ahat0*Xprepre - obj.Bhat0*prepreL );
            obj.delta = Xnow - obj.X_tilde;
        end
        
        function onlineLeastSquare(obj)
            %%% this least square is on the estimation part kind of
            if length(obj.LIP.xf_list) > 3
                actPos = obj.LIP.xf_list(2:end)';
                estPos = obj.MPC_Lresults(1:end,5);
                
                actVel = obj.LIP.dxf_list(2:end)';
                estVel = obj.MPC_Lresults(1:end,6);
                
                deltaPos = actPos - estPos;
                deltaVel = actVel - estVel;
                %%% regress argumenets
                xSSP0 = obj.MPC_Lresults(1:end, 7);
                dxSSP0 = obj.MPC_Lresults(1:end, 8);
                N = length(xSSP0);
                %             paramPos = pinv([xSSP0, dxSSP0, ones(N,1)])*deltaPos;
                %             paramVel = pinv([xSSP0, dxSSP0, ones(N,1)])*deltaVel;
                %             %%% verify
                %             estUpdatePos = estPos' + paramPos'*[xSSP0'; dxSSP0'; ones(1,N)];
                %             estUpdateVel = estVel' + paramVel'*[xSSP0'; dxSSP0'; ones(1,N)];
                %
                paramPos = pinv([xSSP0, dxSSP0])*deltaPos;
                paramVel = pinv([xSSP0, dxSSP0])*deltaVel;
                obj.regreeEst = [paramPos'; paramVel'];
                %%% verify
                estUpdatePos = estPos' + paramPos'*[xSSP0'; dxSSP0'];
                estUpdateVel = estVel' + paramVel'*[xSSP0'; dxSSP0'];
                
                figure(1314159),
                h1 = subplot(2,1,1); hold(h1, 'on');
                plot(1:N, estUpdatePos, 'r', 1:N, actPos, 'b', 1:N, estPos, 'g')
                legend('updateEst', 'act', 'LIP-est')
                title('SSP final position')
                h2 = subplot(2,1,2); hold(h2, 'on');
                plot(1:N, estUpdateVel, 'r', 1:N, actVel, 'b', 1:N, estVel, 'g')
                legend('updateEst', 'act', 'LIP-est')
                title('SSP final velocity'),
            end
        end
        
        function onlineLeastSquareonC(obj, Xi, Li)
            if size(obj.MPC_Lresults, 1) >= 1
                %%
                Xim1 = obj.MPC_Lresults(end, 5:6)';
                Lim1 = obj.MPC_Lresults(end, 1);
                Xi_est = obj.MPC_Lresults(end, end-1:end)';
                %%% Xi_est = Ahat*Xim1 + Bhat*Lim1 + obj.C*Li
                if Li ~=0
                    error = Xi -(Xi_est - obj.C*Li); %% C*Li + delta
                    obj.LIP2SLIP_deltaVec = [obj.LIP2SLIP_deltaVec, error];
                    Lvec = obj.MPC_Lresults(:,1);
                    %%%%[L1, 1; L2, 1; ...; Ln, 1] * [C(1); delta(1)] = obj.LIP2SLIP_deltaVec(1,:);
                    %%%% least square Ax = b
                    A = [Lvec, ones(length(Lvec),1)];
                    b1 = transpose(obj.LIP2SLIP_deltaVec(1, :));
                    paramPos = pinv(A)*b1;
                    
                    %%%%[L1, 1; L2, 1; ...; Ln, 1] * [C(2); delta(2)] = obj.LIP2SLIP_deltaVec(2,:);
                    b2 = transpose(obj.LIP2SLIP_deltaVec(2, :));
                    paramVel = pinv(A)*b2;
                    
                    obj.C = [paramPos(1); paramVel(1)];
                    obj.LIP2SLIP_delta = [paramPos(2); paramVel(2)];
                    
                    %obj.C = (Xi -(Xi_est - obj.C*Li) )/Li;
                else
                    obj.C = zeros(2,1);
                end
            end
        end
        
        function checkParameterConvergence(obj)
            figure,
            for i = 1:3
                for j = 1:2
                    ii = i + (j-1)*3;
                    subplot(2,3,ii), plot( obj.ThetaIDlist(:, ii));
                end
            end
            %% check the gains
            figure, hold on,
            plot(obj.Klist(1,:), 'r')
            plot(obj.Klist(2,:), 'b')
            title('state-feedback gains, deadbeat')
        end
        
        function updateLIPdynamics(obj)
            %%%% update the feedback gains based on the current S2S
            %%%% dynamics
            %%%%
            a = obj.Anow(1,1); b = obj.Anow(1,2); c = obj.Anow(2,1); d = obj.Anow(2,2);
            e = obj.Bnow(1); f = obj.Bnow(2);
            k1 = (a*c*e + c*d*e - a^2*f - b*c*f)/(-c*e^2 + a*e*f - d*e*f + b*f^2);
            k2 = (b*c*e + d^2*e - a*b*f - b*d*f)/(-c*e^2 + a*e*f - d*e*f + b*f^2);
            
            obj.Kdeadbeat = [k1, k2];
            
            k1 = 1/(-1 + a + b*c + d - a*d);
            k2 = -((-c*e*(-a^2*(-1 + d) + d^2 + b*c*(1 + d) + a*(b*c+d-d^2) + e) + (a^2*b*c - a^3*(-1 + d) + b*c*(b*c + d) - d*e + ...
                a*(-b*c*(-2 + d) + e))*f + b* f^2)/((-1 + a + b*c + d - a*d)*(-c *e^2 + f*(a*e - d*e + b*f))));
            k3 = -(((b^2*c^2 + d^3 + b*c*d*(2+d)-a*(b*c*(-1 + d) + d^3))*e-b*(-a^2*(-1 + d) + d^2 + b*c*(1 + d) +a*(b*c+d-d^2))*f)/((-1 + ...
                a + b*c + d - a*d)*(c*e^2 - f*(a*e - d*e + b*f))));
            obj.K3deadbeat = [k1, k2, k3];
            
            obj.A3 =  [1, a-1, b;
                0, a, b;
                0, c, d];
            obj.B3 = [e+1;
                e;
                f];
        end

    end
    
    methods %%% virtual constraint type step size controller,
        %%% should be in closed form or fast to calculate,
        %%% so that this planning is in the control loop
        function [L, dL] = LIPbasedController(obj, t, x, x2f, dx, uNow, varargin)
            obj.HLIP.t = [obj.HLIP.t, t];
            obj.HLIP.x = [obj.HLIP.x, x];
            obj.HLIP.x2f = [obj.HLIP.x2f, x2f];
            obj.HLIP.dx = [obj.HLIP.dx, dx];
            obj.HLIP.u = [obj.HLIP.u, uNow];

            if ~isempty(varargin) %%% if use the estimated states 
                tEst =  varargin{1};  
                [~, x2fpred, dxpred] = obj.LIPbasedEst(x, x2f,dx, tEst);%%%% estimate preimpact state
                %[x, x2f, dx] = obj.LIPbasedEst(x, x2f,dx, tEst);%%%% estimate preimpact state
                obj.HLIP.pred_p = [obj.HLIP.pred_p, x2fpred]; 
                obj.HLIP.pred_v = [obj.HLIP.pred_v, dxpred]; 
            end 
            if length(x2f) == 1
                %%%% instaneout calculation
                if strcmp(obj.whichController, 'Deadbeat')
                    [L, dL] = obj.DeadbeatP1P2(x2f, dx);
                 elseif strcmp(obj.whichController,'DeadbeatParallel')
                    [L, dL] = obj.DeadbeatP1P2withLIPinparallel(x2f, dx);
                elseif strcmp(obj.whichController, 'DiscreteLQR_P1P2')
                    [L, dL] = obj.DiscreteLQR_P1P2(x2f, dx);
                 elseif strcmp( obj.whichController, 'CLF')
                    [L, dL] = obj.LTI_CLF(x2f, dx);
                elseif strcmp( obj.whichController, 'DeadbeatID')
                    [L, dL] = obj.parameterIDwithDeadbeat(x2f, dx);
                elseif strcmp( obj.whichController, 'DiscreteLQR')
                    [L, dL] = obj.DiscreteLQR(x2f, dx);
                elseif strcmp( obj.whichController, 'Robust')
                    [L, dL] = obj.robustLyapunov(x2f, dx);
                elseif strcmp( obj.whichController, 'VSS')
                    [L, dL] = obj.VSS(x2f, dx);
                elseif strcmp( obj.whichController, 'RobustMPC')
                    [L, dL] = obj.robustMPC();
                 elseif strcmp( obj.whichController, 'SLS')
                    [L, dL] = obj.systemLevelSynthesis(x2f, dx);
                else %%% 'StepSequence' %% assume all kinds of MPC
                    [L, dL] = obj.LIPsteppingSequence(x, x2f, dx);
                end

            else %%%%%%% recalculate for state augmentation, should fix later. %%%%%
                if strcmp(obj.whichController, 'Deadbeat') 
                    [L, dL] = obj.DeadbeatP1P2(x2f, dx);
 
                elseif strcmp(obj.whichController,'DeadbeatParallel')
                    [L, dL] = obj.DeadbeatP1P2withLIPinparallel(x2f, dx);
                elseif strcmp(obj.whichController, 'DiscreteLQR_P1P2')
                    [L, dL] = obj.DiscreteLQR_P1P2(x2f, dx);
                elseif strcmp(obj.whichController, 'CLF')
                    [L, dL] = obj.LTI_CLF(x2f(end), dx(end));
                elseif strcmp(obj.whichController, 'DeadbeatID')
                    [L, dL] = obj.parameterIDwithDeadbeat(x2f(end), dx(end));
                elseif strcmp(obj.whichController, 'DiscreteLQR')
                    [L, dL] = obj.DiscreteLQR( x2f(end), dx(end));
                elseif strcmp(obj.whichController, 'Robust')
                    [L, dL] = obj.robustLyapunov( x2f, dx(end));
                elseif strcmp(obj.whichController, 'VSS')
                    [L, dL] = obj.VSS( x2f, dx);
                elseif strcmp(obj.whichController, 'RobustMPC')
                    [L, dL] = obj.robustMPC();
                elseif strcmp( obj.whichController, 'SLS')
                    [L, dL] = obj.systemLevelSynthesis(x2f, dx);
                else %if strcmp(whichController, 'StepSequence')
                    [L, dL] = obj.LIPsteppingSequence( x, x2f, dx);
                end
            end
           %  obj.latestStepL = L; 
        end
        
        function [stepLength, dstepLength] = DeadbeatP1P2(obj, x2f, dx)
            xf_est = x2f + obj.Emeasurement(1); % vector operation
            dxf_est = dx + obj.Emeasurement(2);
            ddxf_est = obj.lambda^2*xf_est;
            xnow = [xf_est; dxf_est];          
            dxnow = [dxf_est; obj.LIPacc];

            theta = obj.theta(end); %atan2(0.8, 3); 
            dstepLength = 0 ; %%% not too much use for aSLIP walking
            switch obj.gaitType 
                case 'P1'
                   % obj.Xdes(1) = obj.vxdes/obj.sigma1 + obj.z*sin(theta);
                    S = transpose((obj.Ahat0 - eye(2))*obj.Xdes)*[1;0]/obj.Bhat0(1);
                     %%%% in the form for feedback gain*error
                    K = obj.Kdeadbeat; 
                    stepLength = K*(xnow - obj.Xdes) - S;
                   % stepLength = stepLength*cos(theta); 
                   % stepLength = K*(xnow - obj.XLIP) + obj.targetStepLength;
                    obj.targetStepLength = K*(obj.XLIP - obj.Xdes) - S;
                case 'P2'
                    if obj.stanceLegNum == -1
                        xdes = obj.XLdes;
                    else
                        xdes = obj.XRdes;
                    end
                    stepLengthDes = 2*xdes(1) + obj.TD*xdes(2);
                     K = obj.Kdeadbeat;
                     stepLength = K*(xnow - xdes) + stepLengthDes;
                    obj.targetStepLength = K*(obj.XLIP - xdes) + stepLengthDes;
                otherwise
                    disp('it should not be here')
            end
            dstepLength = K*dxnow; 
        end 
        
        function [stepLength, dstepLength] = DeadbeatP1P2withLIPinparallel(obj, x2f, dx)
            xf_est =  x2f + obj.Emeasurement(1); % vector operation
            dxf_est = dx + obj.Emeasurement(2);
            ddxf_est = obj.lambda^2*xf_est;
            xnow = [xf_est; dxf_est];
            dxnow = [dxf_est; obj.LIPacc];
            dstepLength = 0 ; %%% not too much use for aSLIP walking
            switch obj.gaitType 
                case 'P1'
                   % obj.Xdes(1) = obj.vxdes/obj.sigma1 + obj.z*sin(theta);
                    S = transpose((obj.Ahat0 - eye(2))*obj.Xdes)*[1;0]/obj.Bhat0(1);
                    % S = vxdes/([0,1]*inv(obj.Anow - eye(2))*obj.Bnow);  %%% the regressed formulation %%%%%%%%%
                    %%%% in the form for feedback gain*error
                    K = obj.Kdeadbeat;  % 1.0000    0.5240
                    stepLength = K*(xnow - obj.XLIP) + obj.targetStepLength;
                    obj.targetStepLength = obj.Kdeadbeat*(obj.XLIP - obj.Xdes) - S;
                case 'P2'
                    if obj.stanceLegNum == -1
                        xdes = obj.XLdes;
                    else
                        xdes = obj.XRdes;
                    end
                    stepLengthDes = 2*xdes(1) + obj.TD*xdes(2);
                    K = obj.Kdeadbeat;
                    stepLength = K*(xnow - xdes) + stepLengthDes;
                    stepLength = K*(xnow - obj.XLIP) + obj.targetStepLength;
                    obj.targetStepLength = obj.Kdeadbeat*(obj.XLIP - xdes) + stepLengthDes;
                otherwise
                    disp('it should not be here')
            end
            dstepLength = K*dxnow; 
        end 
        
        function [stepLength, dstepLength] = robustMPC(obj)
            stepLength = obj.targetStepLength; 
            dstepLength = 0*stepLength;
        end
        
        function [stepLength, dstepLength] = LIPsteppingSequence(obj, x, x2f, dx)
            xnow = [x2f; dx];
            x3now = [x; xnow];
            
            K = obj.Kdeadbeat; 
             switch obj.goalBehavior 
                case 'local'
                    stepLength = obj.targetStepLength + K*(xnow - obj.XLIP); %% deadbeat gain
                case 'global'
                    stepLength = obj.targetStepLength + K*(x3now - obj.X3LIP); %% deadbeat gain
                   % stepLength = obj.targetStepLength + obj.K3LQR*(x3now - obj.X3LIP);
                otherwise
                    disp('wrong behavior')
             end 
             dstepLength = 0*stepLength;
        end
        
        function [stepLength, dstepLength] = parameterIDwithDeadbeat(obj, varargin)
            if isempty(varargin)
                xf_est = obj.polar.x2f;
                dxf_est = obj.polar.dx;
            else
                xf_est = varargin{1}; % vector operation
                dxf_est = double( varargin{2});
            end
            xdes = obj.vxdes/obj.sigma1;
            Xdes = [xdes; obj.vxdes];
            
            %%%% get the updated dynamics from parameter ID
            thetaNow = transpose(obj.ThetaID);
            Anow = thetaNow(:, 1:2) + obj.Am;
            Bnow = thetaNow(:, 3);
            
            %%%% shift the desired state based on the identified A, B %%%%%
            temp = (Anow - eye(2))\Bnow;
            xdes = obj.vxdes/temp(2)*temp(1);
            Xdes = [xdes; obj.vxdes];
            
            %%% state feedback, deadbeat controller gains.
            a = Anow(1,1); b = Anow(1,2); c = Anow(2,1); d = Anow(2,2);
            e = Bnow(1); f = Bnow(2);
            k1 = (a*c*e + c*d*e - a^2*f - b*c*f)/(-c*e^2 + a*e*f - d*e*f + b*f^2);
            k2 = (b*c*e + d^2*e - a*b*f - b*d*f)/(-c*e^2 + a*e*f - d*e*f + b*f^2);
            
            L = k1*(xf_est - Xdes(1)) + k2*(dxf_est - Xdes(2));
            
            S = transpose((Anow -eye(2))*Xdes)*[1;0]/Bnow(1);
            stepLength = L - S;
            
            dstepLength = 0*dxf_est;
          
            obj.Know = [k1; k2];
        end
        
        function [stepLength, dstepLength] = systemLevelSynthesis(obj, x2f, dx)
            xf_est = x2f + obj.Emeasurement(1); % vector operation
            dxf_est = dx + obj.Emeasurement(2);
            ddxf_est = obj.lambda^2*xf_est;
            xnow = [xf_est; dxf_est];
            dxnow = [dxf_est; obj.LIPacc];
            
           %  dstepLength = 0; %%% not too much use for aSLIP walking
            switch obj.gaitType
                case 'P1'
                    [Usls, dUsls] = obj.calculateSLS(xnow - obj.XLIP, dxnow);
                    % S = transpose((obj.Ahat0 - eye(2))*obj.Xdes)*[1;0]/obj.Bhat0(1);
                    stepLength = Usls + obj.targetStepLength;
                    obj.targetStepLength = obj.Kdeadbeat*(obj.XLIP - obj.Xdes) + obj.Udes;
                    obj.SLS_ulog = [obj.SLS_ulog, Usls];
                case 'P2'
                    if obj.stanceLegNum == -1
                        xdes = obj.XLdes;
                    else
                        xdes = obj.XRdes;
                    end
                    stepLengthDes = 2*xdes(1) + obj.TD*xdes(2);
                    obj.targetStepLength = stepLengthDes;
                    Usls = obj.calculateSLS(xnow - xdes);
                    stepLength = Usls + stepLengthDes;
            end
            dstepLength = dUsls; 
        end
        
        function getSLSgain(obj)
            A = obj.Ahat0;
            B = obj.Bhat0;
            C = [eye(2)*50; zeros(1,2)];
            D = [zeros(2,1); 1];
            %%%%%%%%%%%%%% |C x + D u |^2
            Nx = 2;
            Nu = 1;
            T = obj.N_fir; %%%%% FIR steps to zero
            
            Nw = 2;
            cvx_begin
            cvx_precision high
            
            % cvx_solver sedumi
            % expression Delta(Nx,Nx*T)
            variable Rs(Nx,Nw,T) %%% pfi_x
            variable Ms(Nu,Nw,T) %%% pfi_u
            
            %populate decision variables
            %locality constraints automatically enforced by limiting support of R and M
            for t = 1:T
                R{t} = Rs(:,:,t);
                M{t} = Ms(:,:,t);
            end
            lambda = 1;
            objective = compute_H2(R,M,C,D,T);
            R{1} == eye(Nx);
            
            for t= 1:T-1
                %         Delta(:,(t-1)*Nx+1:t*Nx) = R{t+1} - A*R{t} - B*M{t};
                R{t+1} == A*R{t} + B*M{t};
            end
            R{T} == zeros(Nx,Nx);
            
            %     robust_stab = norm(Delta,inf);
            
            % minimize(objective+lambda*robust_stab)
            minimize(objective)
            cvx_end
            
            obj.SLS_Ugain = M; 
            obj.SLS_XrefGain = R; 
        end
        
        function getSLSgainwithConst(obj)
            A = obj.Ahat0;
            B = obj.Bhat0;
            C = [eye(2); zeros(1,2)];
            D = [zeros(2,1); 10];
            %%%%%%%%%%%%%% |C x + D u |^2
            ep_max = 0.2; ev_max = 1.3; %%% external disturbance
            dp_max = 0.1; dv_max = 0.1; 
            p_max = 0.8;
            p_min = -.8;
            v_bound = 2; %%% can be larger range
            ubound = 0.4; 
            u_min = -ubound - obj.Udes;
            u_max = ubound - obj.Udes; 
            d = zeros(2,1); 
            Nx = 2;
            Nu = 1;
            T = obj.N_fir; %%%%% FIR steps to zero
            
            Nw = 2;
            
            %%% configure constraints
            Hx = kron(eye(Nx), [1;-1]);  %%%% Hx < h
            Hu = kron(eye(Nu), [1;-1]);  %%%% Hu < hu
            G = kron(eye(Nx), [1;-1]);
         
            h = [p_max; -p_min; v_bound; v_bound];
            hu = [u_max; -u_min];
            
            g_vec{1} = [ep_max; ep_max; ev_max; dv_max];
            for t = 2:T
                g_vec{t} = [dp_max; dp_max; dv_max; dv_max];
            end

            [px,~] = size(Hx);
            [pu,~] = size(Hu);
            [q,~] = size(G);

            cvx_begin
            cvx_precision high
            
            % cvx_solver sedumi
            % expression Delta(Nx,Nx*T)
            variable Rs(Nx,Nw,T) %%% pfi_x
            variable Ms(Nu,Nw,T) %%% pfi_u
            variable lambdas_x(px,q,T*(T+1)/2)
            variable lambdas_u(pu,q,T*(T+1)/2)
            %populate decision variables
            %locality constraints automatically enforced by limiting support of R and M
            lambda_count = 1;
            for t = 1:T
                R{t} = Rs(:,:,t);
                M{t} = Ms(:,:,t);
                for t1=1:t
                    lambda_x{t,t1} = lambdas_x(:,:,lambda_count);
                    lambda_u{t,t1} = lambdas_u(:,:,lambda_count);
                    lambda_count = lambda_count+1;
                end
            end
            objective = compute_H2(R,M,C,D,T);
            R{1} == eye(Nx);
            
            for t= 1:T-1
                %         Delta(:,(t-1)*Nx+1:t*Nx) = R{t+1} - A*R{t} - B*M{t};
                R{t+1} == A*R{t} + B*M{t} + d;
            end
            R{T} == zeros(Nx,Nx);
            for t=1:T
                lambda_vec=[];
                lambda_gx=0*h;
                lambda_gu=0*hu;
                for t1=1:t
                    Hx*R{t1} == lambda_x{t,t1}*G;
                    lambda_gx = lambda_gx+lambda_x{t,t1}*g_vec{T-t1+1};
                    lambda_x{t,t1}>=0
                    
                    Hu*M{t1} == lambda_u{t,t1}*G;
                    lambda_gu = lambda_gu+lambda_u{t,t1}*g_vec{T-t1+1};
                    lambda_u{t,t1}>=0
                end
                
                lambda_gx<=h;
                lambda_gu<=hu;
            end
            %     robust_stab = norm(Delta,inf);
            
            % minimize(objective+lambda*robust_stab)
            minimize(objective)
            cvx_end
            
            obj.SLS_Ugain = M;
            obj.SLS_XrefGain = R;
        end
        
        function getSLSgainwithRS2S(obj, Fpush, uMax)
            A = obj.Abar;
            B = obj.Bbar;
            Ccost = [eye(2); zeros(1,2)];
            Dcost = [zeros(2,1); 1];
            %%%%%%%%%%%%%% |C x + D u |^2
            ep_max = 0.18; ev_max = 1.2; %%% external disturbance (to do: derive it again)
            obj.getPushDisturbance(Fpush);
            ep_max = obj.pushDisturbance(1);
            ev_max = obj.pushDisturbance(2);
            dp_max = obj.episonMax(1); dv_max = obj.episonMax(2); %% model
            p_max = 1;
            p_min = -1;
            v_bound = 10; %%% can be larger range
            ubound = uMax;%  0.75;
            u_min = -ubound - obj.Udes;
            u_max = ubound - obj.Udes;
          %   ubound = identifyUmaxFromPush(obj);
            Nx = 2;
            Nu = 1; Nw = 2; 
            T = obj.N_fir; %%%%% FIR steps to zero
            %%% configure constraints
            Hx = kron(eye(Nx), [1;-1]);  %%%% Hx < h
            Hu = kron(eye(Nu), [1;-1]);  %%%% Hu < hu
            G = kron(eye(Nx), [1;-1]);
            
            h = [p_max; -p_min; v_bound; v_bound];  %%% state constraint
            hu = [u_max; -u_min];
            NN = 3; 
            g_vec{1} = NN*[dp_max; dp_max; dv_max; dv_max];
            g_vec{2} = [ep_max+dp_max; ep_max+dp_max; ev_max+dv_max;ev_max+dv_max];
            for t = 3:T
                g_vec{t} = [dp_max; dp_max; dv_max; dv_max];
            end
            [px,~] = size(Hx);
            [pu,~] = size(Hu);
            [q,~] = size(G);
            cvx_begin
            cvx_precision high
            % cvx_solver sedumi
            % expression Delta(Nx,Nx*T)
            variable Rs(Nx,Nw,T) %%% pfi_x
            variable Ms(Nu,Nw,T) %%% pfi_u
            variable lambdas_x(px,q,T*(T+1)/2)
            variable lambdas_u(pu,q,T*(T+1)/2)
            %populate decision variables %locality constraints automatically enforced by limiting support of R and M
            lambda_count = 1;
            for t = 1:T
                R{t} = Rs(:,:,t);
                M{t} = Ms(:,:,t);
                for t1=1:t
                    lambda_x{t,t1} = lambdas_x(:,:,lambda_count);
                    lambda_u{t,t1} = lambdas_u(:,:,lambda_count);
                    lambda_count = lambda_count+1;
                end
            end
            objective = compute_H2(R,M,Ccost,Dcost,T);
            R{1} == eye(Nx);
            
            for t= 1:T-1
                %         Delta(:,(t-1)*Nx+1:t*Nx) = R{t+1} - A*R{t} - B*M{t};
                R{t+1} == A*R{t} + B*M{t};
            end
            R{T} == zeros(Nx,Nx);
            for t=1:T
                lambda_vec=[];
                lambda_gx=0*h;
                lambda_gu=0*hu;
                for t1=1:t
                    Hx*R{t1} == lambda_x{t,t1}*G;
                    lambda_gx = lambda_gx+lambda_x{t,t1}*g_vec{t-t1+1};
                    lambda_x{t,t1}>=0
                    
                    Hu*M{t1} == lambda_u{t,t1}*G;
                    lambda_gu = lambda_gu+lambda_u{t,t1}*g_vec{t-t1+1};
                    lambda_u{t,t1}>=0
                end
                if  t <  T
                    lambda_gx<=h;
                else
                    lambda_gx <= NN*[dp_max;dp_max;dv_max; dv_max];
                end
               lambda_gu<=hu;
            end
            %     robust_stab = norm(Delta,inf);
           %  minimize(objective+lambda*robust_stab)
           
            minimize(objective)
            cvx_end
            
            obj.SLS_Ugain = M;
            obj.SLS_XrefGain = R;
        end
       
        function uMax = identifyFromPush(obj, Fpush) 
            A = obj.Abar;
            B = obj.Bbar;
            obj.getPushDisturbance(Fpush);
            ep_max = obj.pushDisturbance(1);
            ev_max = obj.pushDisturbance(2);
            dp_max = obj.episonMax(1); dv_max = obj.episonMax(2); %% model
            p_max = 0.6;
            p_min = -0.6;
            v_bound = 3; %%% can be larger range
            Nx = 2;
            Nu = 1;
            T = obj.N_fir; %%%%% FIR steps to zero
            Nw = 2;
            %%% configure constraints
            Hx = kron(eye(Nx), [1;-1]);  %%%% Hx < h
            Hu = kron(eye(Nu), [1;-1]);  %%%% Hu < hu
            G = kron(eye(Nx), [1;-1]);
           
            h = [p_max; -p_min; v_bound; v_bound];  %%% state constraint
            
            g_vec{1} = [ep_max+dp_max; ep_max+dp_max; ev_max+dv_max;ev_max+dv_max];
            for t = 2:T
                g_vec{t} = [dp_max; dp_max; dv_max; dv_max];
            end
            [px,~] = size(Hx);
            [pu,~] = size(Hu);
            [q,~] = size(G);
            cvx_begin
            cvx_precision high 
            variable uMax(1) %%% pfi_x
            u_min = -uMax - obj.Udes;
            u_max = uMax - obj.Udes;
            hu = [u_max; -u_min];

            % cvx_solver sedumi
            % expression Delta(Nx,Nx*T)
            variable Rs(Nx,Nw,T) %%% pfi_x
            variable wMax(Nw) %%% pfi_x
            variable Ms(Nu,Nw,T) %%% pfi_u
            variable lambdas_x(px,q,T*(T+1)/2)
            variable lambdas_u(pu,q,T*(T+1)/2)
            %populate decision variables %locality constraints automatically enforced by limiting support of R and M
            lambda_count = 1;
            for t = 1:T
                R{t} = Rs(:,:,t);
                M{t} = Ms(:,:,t);
                for t1=1:t
                    lambda_x{t,t1} = lambdas_x(:,:,lambda_count);
                    lambda_u{t,t1} = lambdas_u(:,:,lambda_count);
                    lambda_count = lambda_count+1;
                end
            end
           %  objective = compute_H2(R,M,Ccost,Dcost,T);
            objective = uMax;
            R{1} == eye(Nx);
            
            for t= 1:T-1
                %         Delta(:,(t-1)*Nx+1:t*Nx) = R{t+1} - A*R{t} - B*M{t};
                R{t+1} == A*R{t} + B*M{t};
            end
            R{T} == zeros(Nx,Nx);
            for t=1:T
                lambda_vec=[];
                lambda_gx=0*h;
                lambda_gu=0*hu;
                for t1=1:t
                    Hx*R{t1} == lambda_x{t,t1}*G;
                    lambda_gx = lambda_gx+lambda_x{t,t1}*g_vec{t-t1+1};
                    lambda_x{t,t1}>=0
                    
                    Hu*M{t1} == lambda_u{t,t1}*G;
                    lambda_gu = lambda_gu+lambda_u{t,t1}*g_vec{t-t1+1};
                    lambda_u{t,t1}>=0
                end
                lambda_gx<=h;
                lambda_gu<=hu;
            end
            %     robust_stab = norm(Delta,inf);
            
            % minimize(objective+lambda*robust_stab)
            minimize(objective)
            cvx_end
          
        end 
        
        function [Usls, dUsls] = calculateSLS(obj, x, dx)
            Usls = 0;
            dUsls = 0; 
            obj.x_ref_k = zeros(2,1); 
            %%%% get previous estimated w (not the real w)
            w_est = x - obj.x_ref_km1;%%% previous  
            obj.w_est_vec(:,end) = w_est; 
            for jj = 1:obj.N_fir
                Usls = Usls + obj.SLS_Ugain{jj}*obj.w_est_vec(:,obj.N_fir+1-jj);   %%%% get the current desired U.
            end
            dUsls = dUsls +  obj.SLS_Ugain{1}*dx;
            % x(:,i+1) = A*x(:,i) + B1*w_d(:,i)+ B*u(:,i); for simulation purpose
            for jj = 2:obj.N_fir
                obj.x_ref_k = obj.x_ref_k + obj.SLS_XrefGain{jj}*obj.w_est_vec(:,obj.N_fir+2-jj);
            end
        end 
        
        function ubound = identifyUmaxFromPush(obj) 
            A = obj.Abar;
            B = obj.Bbar;
            ep_max = obj.pushDisturbance(1);
            ev_max = obj.pushDisturbance(2);
            dp_max = obj.episonMax(1); dv_max = obj.episonMax(2); %% model
            p_max = 0.5;  p_min = -0.5;
            v_bound = 10; %%% can be larger range
            
            N = obj.N_fir;
            ubound = 1;
            disb = sdpvar(1); 
            x = sdpvar(2, N);
            u = sdpvar(1, N); 
            w = sdpvar(2,N); 
            u_min = -ubound - obj.Udes;
            u_max = ubound - obj.Udes;
            
            %%% argmin ubound
            %%%% u 
            F = []; C = []; W = []; 
            W = [W, x(:,1) == zeros(2,1)];
            W = [W, x(:,N) == zeros(2,1)];

            W = [W, - disb*0.3<= w(1,1)<= disb*0.3];
            W = [W, - disb <= w(2,1) <= disb];
            for i = 1:N-1 
                F = [F,  x(:,i+1) == A*x(:,i) + B*u(i) + w(:,i)];
                W = [W,  -dp_max <= w(1,1+i) <= dp_max];
                W = [W,  -dv_max <= w(2,1+i) <= dv_max];
            end
            for i = 1:N
                F = [F, u_min <= u(i) <= u_max];
                W = [W, p_min <= x(1, i) <= p_max];
                W = [W, -v_bound <= x(2, i) <= v_bound];
            end
           %  W = [W, uncertain(w)];

            objective = - disb^2 ; 
            optimize([F, W],objective);
            disb = value(disb)
            x = value(x); 
            w = value(w);
            u = value(u)
            
            %%
            ubound = 0.5;
            u = sdpvar(1); 
            w = sdpvar(2,1); 
            wmax = sdpvar(1); 
            F = []; W = [];
            F = [F,   B*u + w ==0];            W = [W, -ubound <= u <= ubound];
        %      F = [F,   B*u + w >= -0.1];

            W = [W, -wmax <=  w(1) <= wmax];
            W = [W, -3*wmax <=  w(2) <= 3*wmax];
           %  W = [W,  uncertain(w)];
            objective = -wmax^2; 
            solvesdp([F, W], objective);
            wmax = value(wmax)
        end 
        
    end
    
    methods %%%% MPC type planning, at discrete level 

        function updateLIPstates(obj, t0)
            %%%%% used for the MPC part, presumably 
            %%%%% the nominal model of reference. 
            ifupdate = 0; %%%%
           
            %%%% SLS update
            obj.x_ref_km1 = obj.x_ref_k; 
            obj.w_est_vec(:,1:end-1) = obj.w_est_vec(:,2:end);
            %%% save robot final state
            if ~isempty(obj.HLIP.t)
                obj.robotXF.t = [obj.robotXF.t, obj.HLIP.t(end)];
                obj.robotXF.p = [obj.robotXF.p, obj.HLIP.x2f(end)];
                obj.robotXF.x = [obj.robotXF.x, obj.HLIP.x(end)];
                obj.robotXF.v = [obj.robotXF.v, obj.HLIP.dx(end)];
                obj.robotXF.u = [obj.robotXF.u, obj.HLIP.u(end)];
                obj.XLIPlist = [obj.XLIPlist, obj.XLIP];
                obj.X3LIPlist = [obj.X3LIPlist, obj.X3LIP];
                obj.targetStepLengthVec = [obj.targetStepLengthVec, obj.targetStepLength]; 
            end
            
            if ifupdate == 1 %%% if update the linear dynamics
                obj.updateLIPdynamics();    %%%% if want to update the deadbeat gain
                obj.XLIP = obj.Anow*obj.XLIP + obj.Bnow*obj.targetStepLength;
                obj.X3LIP = obj.A3*obj.X3LIP + obj.B3*obj.targetStepLength;
                %%%% if not use the parameter ID
            elseif obj.ifRS2S == 1
                obj.XLIP = obj.Abar*obj.XLIP + obj.Bbar*obj.targetStepLength+ obj.Cbar;
            else 
                obj.XLIP = obj.Ahat0*obj.XLIP + obj.Bhat0*obj.targetStepLength;
                obj.X3LIP = obj.A3*obj.X3LIP + obj.B3*obj.targetStepLength;
            end
        end
        
        function MPC_controller(obj, type, t0, Xnow) 
                % obj.LMPC_LIP();
                switch type
                    case 'MPCtrackingVelocityP1'
                         obj.MPC_P1();
                    case 'MPCtrackingVelocityP2'
                         obj.MPC_P2();
                    case 'MPCtrackingPath'
                        obj.MPC_LIP_path(t0);
                    case 'MPCwithBarrier'
                        obj.MPC_LIP_barrier(t0);
                    case 'RobustMPC'
                        obj.MPC_Robust(t0, Xnow); 
                end
        end 
        
        function MPC_Robust(obj, t, X3now)
            %%% test for fixed behavior %%%%
            Nsteps = 6; %%% 4 or 6 or 10 is the same, smaller is not good..
            switch obj.goalBehavior
                case 'local'
                    Xnow = X3now(2:3); 
                    K = obj.KposDeadbeat; 
                    Q = eye(2)*1e+2; R = 1;
                    [ x_cl, stepvec ] = solve_CFTOCP( Xnow, obj.Xdes, Nsteps, Q, R, obj.Ahat0, obj.Bhat0, ...
                        obj.Xset, obj.Uset, 'quadprog', obj.Eset, K);
                    obj.XLIP = x_cl(:, 1);                    
                case 'global'             %%% 3 states tracking path
                    Q = blkdiag(1, 0, 1); R = 1;
                    K = obj.K3deadbeat;
                    Xnow = X3now; 
                    [ x_cl, stepvec ] = obj.optStepping.robustMPC(t, obj.TS+obj.TD, Xnow, Nsteps, Q, R, obj.A3, obj.B3, ...
                        obj.X3set, obj.Uset, obj.E3set, K);
                    obj.X3LIP = x_cl(:, 1);
                otherwise
                    disp('wrong behavior');
            end
            
           % obj.targetStepLength = stepvec(1);
             obj.targetStepLength = stepvec(1) + K*(Xnow - x_cl(:,1));
            obj.MPCsol.usol = [obj.MPCsol.usol; [t, stepvec] ]; 
            [nrow, ncol] = size(x_cl);
            xsol = zeros(nrow, ncol, 2); 
            xsol(:, :, 1) = x_cl; 
            obj.MPCsol.xsol = [obj.MPCsol.xsol; xsol(:,:,1)]; 
        end
        
        function MPC_P1(obj)
            
            Nsteps = 7;
            vxdes0 = obj.vxdes; %obj.massTraj.SSP.dx(end);
            DistSum = vxdes0*(obj.TS + obj.TD); % total traveled distance in SSP and DSP.
            vxdes = DistSum/(2/obj.sigma1 + obj.TD);
            
            obj.Xdes = [vxdes/obj.sigma1;
                        vxdes];
            
            [stepVec, x_cl] =  obj.optStepping.MPConVelocityConciseP1(obj.TD, obj.TS, Nsteps, obj.XLIP, obj.Ahat0, obj.Bhat0, obj.Xdes);
            obj.targetStepLength = stepVec(1);
            
           % obj.XLIP = obj.Ahat0*obj.XLIP + obj.Bhat0*stepVec(1);
        end
        
        function MPC_P2(obj)
            
            Nsteps = 4;
            d = obj.lambda^2*(sech(obj.lambda*obj.TS/2))^2*obj.vxdes*(obj.TS+obj.TD)/(obj.lambda^2*obj.TD + 2*obj.sigma2);
            X1Fdes = -obj.targetStepSize/(2 + obj.TD*obj.sigma2);
            obj.x1Fdes = X1Fdes;
            VX1Fdes = obj.sigma2*X1Fdes + d;
            obj.XLdes = [obj.x1Fdes; VX1Fdes];
            [x2Fdes, VX2Fdes, ~] = LIP_sol(obj.lambda, -X1Fdes, VX1Fdes, obj.TS);
            obj.XRdes = [x2Fdes; VX2Fdes];
            
            [stepVec, x_cl] =  obj.optStepping.MPConVelocityConciseP2(obj.TD, obj.TS, Nsteps, obj.XLdes,...
                obj.Ahat0, obj.Bhat0,  obj.stanceLegNum, obj.XLdes, obj.XRdes);
            obj.targetStepLength = stepVec(1);
            
            % obj.XLIP = obj.Ahat0*obj.XLIP + obj.Bhat0*stepVec(1);
        end
        
        function MPC_LIP(obj)
            Lnorm0 = obj.stepLengthSequence(end);
            
            % if called at the end of SSP, planning for the next step length
            Xprev = [obj.polar.x2f; obj.polar.dx]; % this is the state at SSP- or DSP+, the current step size is already fix.
            % estimate the next step SSP- state from the DSP+
            xSSP0 = - Lnorm0 + Xprev(1) + obj.TD*Xprev(2);
            dxSSP0 = Xprev(2);
            
             [xSSPF, dxSSPF] = obj.LIPsolution(xSSP0, dxSSP0, obj.TS);
            X0F = (eye(2))*[xSSPF; dxSSPF];
 
%              X0F = (eye(2)+regreeEst)*[xSSPF; dxSSPF];

             % obj.LIP.est_LIP_states = X0F; 
             %%%%%%
             Nsteps = 4; 
             xdes = obj.vxdes/obj.sigma1;
             Xdes = [xdes; obj.vxdes]; 
             LnormF = xdes + obj.TD*obj.vxdes + obj.vxdes/obj.sigma1;
              X0F = obj.XLIP; 
              Xprev = inv(obj.Anow)*(obj.XLIP - obj.Bnow*obj.targetStepLength); 

            [stepVec, X1LIP] = obj.optStepping.includeStateVariable(obj.TD, obj.TS, obj.lambda, Nsteps, X0F, Xdes, obj.LIP2SLIP_delta, obj.Anow, obj.Bnow);
 
             obj.LIP.est_LIP_states = [obj.LIP.est_LIP_states, X1LIP]; 
             
             if isempty(stepVec)
                 obj.targetStepLength = Lnorm0;
                 addpendthis = [ zeros(1,Nsteps), X0F', xSSP0, dxSSP0, X1LIP' ];
             else
                 obj.targetStepLength = stepVec(1);
                 addpendthis = [stepVec', X0F',   xSSP0, dxSSP0, X1LIP'];
             end
             obj.MPC_Lresults = [obj.MPC_Lresults; addpendthis];
             
        end 
        
        function MPC_LIP_path(obj, t)
            Lnorm0 = obj.stepLengthSequence(end);
            % if called at the end of SSP, planning for the next step length
            X0F = [obj.polar.x2f; obj.polar.dx] + obj.Emeasurement; % this is the state at SSP- or DSP+, the current step size is already fix.
            % estimate the next step SSP- state from the DSP+
            xSSP0 = - Lnorm0 + X0F(1) + obj.TD*X0F(2);
            dxSSP0 = X0F(2);
            
            %              %% estimate the SSP- state from SSP+
            %              xSSP0 = obj.LIP.x0actVec(end);
            %              dxSSP0 = obj.polar.dx;
            [xSSPF, dxSSPF] = obj.LIPsolution(xSSP0, dxSSP0, obj.TS);
            X0F = [obj.polar.x-xSSP0+xSSPF + obj.TD*X0F(2); xSSPF; dxSSPF];
            
            %%%%%%
            Nsteps = 6;
            xdes = obj.vxdes/obj.sigma1;
            Xdes = [xdes; obj.vxdes];
            LnormF = xdes + obj.TD*obj.vxdes + obj.vxdes/obj.sigma1;
            
            X0F = obj.X3LIP; 
            [stepVec, X1LIP] = obj.optStepping.MPConPathSagittal(t, obj.TD, obj.TS, obj.lambda, Nsteps, X0F, obj.Anow, obj.Bnow);
            % [stepVec, X1LIP] = obj.optStepping.MPConPathLateral(t, obj.TD, obj.TS, obj.lambda, Nsteps, X0F, length(obj.stepLengthSequence), obj.Anow, obj.Bnow);
            obj.LIP.est_LIP_states = [obj.LIP.est_LIP_states, X1LIP(2:3)];
            
            if isempty(stepVec)
                obj.targetStepLength = Lnorm0;
                addpendthis = [ zeros(1,Nsteps), X0F', xSSP0, dxSSP0, X1LIP' ];
            else
                obj.targetStepLength = stepVec(1);
                addpendthis = [stepVec', X0F',   xSSP0, dxSSP0, X1LIP'];
            end
            obj.MPC_Lresults = [obj.MPC_Lresults; addpendthis];
        end
        
        function MPC_LIP_barrier(obj,t)
            Lnorm0 = obj.stepLengthSequence(end);
            
            % if called at the end of SSP, planning for the next step length
            X0F = [obj.polar.x2f; obj.polar.dx]; % this is the state at SSP- or DSP+, the current step size is already fix.
            % estimate the next step SSP- state from the DSP+
            xSSP0 = - Lnorm0 + X0F(1) + obj.TD*X0F(2);
            dxSSP0 = X0F(2);

            [xSSPF, dxSSPF] = obj.LIPsolution(xSSP0, dxSSP0, obj.TS);
            X0F = [obj.polar.x-xSSP0 + xSSPF + obj.TD*X0F(2); xSSPF; dxSSPF];
            
            %%%%%%
            Nsteps = 4;
            xdes = obj.vxdes/obj.sigma1;
            Xdes = [xdes; obj.vxdes];
            LnormF = xdes + obj.TD*obj.vxdes + obj.vxdes/obj.sigma1;
            [stepVec, X1LIP] = obj.optStepping.MPConPath(t, obj.TD, obj.TS, obj.lambda, Nsteps, X0F);
            obj.LIP.est_LIP_states = [obj.LIP.est_LIP_states, X1LIP(2:3)];
            
            if isempty(stepVec)
                obj.targetStepLength = Lnorm0;
                addpendthis = [ zeros(1,Nsteps), X0F', xSSP0, dxSSP0, X1LIP' ];
            else
                obj.targetStepLength = stepVec(1);
                addpendthis = [stepVec', X0F',   xSSP0, dxSSP0, X1LIP'];
            end
            obj.MPC_Lresults = [obj.MPC_Lresults; addpendthis];
            
        end
        
        function LMPC_LIP(obj)
            %%%% test Lyapunov constraint in the MPC 
            % if called at the end of SSP, planning for the next step length 
             Lnorm0 = obj.stepLengthSequence(end);
             X0F = [obj.polar.x2f; obj.polar.dx]; % this is the state at SSP- or DSP+, the current step size is already fix. 
            % estimate the next step SSP- state from the DSP+ 
             xSSP0 = - Lnorm0 + X0F(1) + obj.TD*X0F(2); 
             dxSSP0 = X0F(2); 
             
%              %% estimate the SSP- state from SSP+
%              xSSP0 = obj.polar.x2f;
%              dxSSP0 = obj.polar.dx;
%              
             [xSSPF, dxSSPF] = obj.LIPsolution(xSSP0, dxSSP0, obj.TS);
             X0F = [xSSPF; dxSSPF];
             % obj.LIP.est_LIP_states = X0F; 
             %%%%%%
             Nsteps = 4; 
             obj.sigma1 = obj.lambda*coth(obj.TS/2*obj.lambda);
             xdes = obj.vxdes/obj.sigma1;
             Xdes = [xdes; obj.vxdes]; 
             [stepVec, X1LIP] = obj.optStepping.MPConError(obj.TD, obj.TS, obj.lambda, Nsteps, X0F, Xdes, obj.LIP2SLIP_delta);
            % [stepVec, X1LIP] = obj.optStepping.MPCwithLaypunovQCQP(TD, TS, obj.lambda, Nsteps, X0F, Xdes, obj.LIP2SLIP_delta);
             obj.LIP.est_LIP_states = [obj.LIP.est_LIP_states, X1LIP]; 
             
             if isempty(stepVec)
                 obj.targetStepLength =  obj.stepLengthSequence(end);
                 addpendthis = [ zeros(1,Nsteps), X0F', xSSP0, dxSSP0, X1LIP' ];
             else
                 obj.targetStepLength = stepVec(1);
                 addpendthis = [stepVec', X0F',   xSSP0, dxSSP0, X1LIP'];
             end
             obj.MPC_Lresults = [obj.MPC_Lresults; addpendthis];
             
        end 
        
        function MPC3states_LIP(obj)
            if ~isempty(obj.stepLengthSequence) 
                Lnorm0 = obj.stepLengthSequence(end);
            else 
                Lnorm0 = 0; 
            end
            Lnorm0 = obj.curstepLength; 
            %%%% 
             X0F = [obj.polar.x2f; obj.polar.dx];
             %%%%%%
             obj.onlineLeastSquareonC(X0F, Lnorm0)
             
             Nsteps = 4; 
             xdes = obj.vxdes/obj.sigma1;
             Xdes = [xdes; obj.vxdes]; 
             LnormF = xdes + obj.TD*obj.vxdes + obj.vxdes/obj.sigma1;
             [stepVec, X1LIP] = obj.optStepping.threeStateFormulation(obj.TD, obj.TS, obj.lambda, Nsteps, X0F, Lnorm0, Xdes, obj.C, obj.LIP2SLIP_delta); 
             obj.LIP.est_LIP_states = [obj.LIP.est_LIP_states, X1LIP]; 
             
             if isempty(stepVec)
                 obj.targetStepLength = Lnorm0;
                 addpendthis = [ zeros(1,Nsteps), X0F', 0, 0, X1LIP' ];
             else
                 obj.targetStepLength = stepVec(1);
                 addpendthis = [stepVec', X0F',   0, 0, X1LIP'];
             end
             obj.MPC_Lresults = [obj.MPC_Lresults; addpendthis];
             
        end 
       
    end 
    
    methods   % support function  
        
        function [xsol, dxsol] = LIPsolution(obj, x0, dx0, t)
            v = obj.lambda;
            c1 = @(y0, yd0) 1/2*(y0 + yd0/v);
            c2 = @(y0, yd0) 1/2*(y0 - yd0/v);
            xsol  =  c1(x0, dx0)*exp(v*t) + c2(x0, dx0)*exp(-v*t);
            dxsol = v*(c1(x0, dx0)*exp(v*t) - c2(x0, dx0)*exp(-v*t));
        end
        
        function calculateNominalTrajPerStep(obj, t0, TS, TD, polarX, curstepLength) 
            %%% given the current state and step size, calculate the HLIP
            %%% trajectory to the next step, 
            %%% use it at preimpact state
            tDSPvec = t0:0.009:t0+TD; 
            tSSPvec = t0+TD:0.009:t0+TD+TS; 
            
            dxDSPvec = ones(size(tDSPvec))*polarX.dx; 
            xDSPvec = polarX.x + cumsum(dxDSPvec)*0.009; 
            
            pSSP0 = - curstepLength + polarX.x2f + polarX.dx*obj.TD; 
            vSSP0 = polarX.dx; 
            [xSSPvec, dxSSPvec] = obj.LIPsolution(pSSP0, vSSP0, tSSPvec-t0-TD);
            xSSPvec = xSSPvec + polarX.x - polarX.x2f + curstepLength; 
            
            obj.NominalTraj.t = [obj.NominalTraj.t, tDSPvec, tSSPvec]; 
            obj.NominalTraj.x = [obj.NominalTraj.x, xDSPvec, xSSPvec]; 
            obj.NominalTraj.dx = [obj.NominalTraj.dx, dxDSPvec, dxSSPvec]; 
        end 
 
        function [ invariantSet ] = computeInvariant(obj, Acl, W)
             %computeInvariant: If the algorithm in this function converges it returns a Robust Positive Invariant (RPI) set
            %   Acl: This matrix is the closed-loop matrix given by Acl = A - BK
            %   W: This polytope defined the domain of the uncertanty. It is defined
            %   using MPC. An 2d example is W = Polyhedron([1 1;1 -1; -1 1; -1 -1]);
            
            X{1} = zeros(size(Acl,1),1); % Initialize the set
            for i = 1:10000
                X{i+1} = Acl*X{i} + W; % Propagate the uncertanty
                X{i+1}.minVRep() % Compute minimal representation
                X{i+1}.minHRep() % Compute minimal representation
                
                % Check if the algorithm has covnerged
                if i > 1
                    if (X{i+1}.contains(X{i})) && (X{i}.contains(X{i+1}))
                        invariantSet = X{i+1}; % Set invaraint to the current iterate
                        disp(['Invariant set computed in i = ',num2str(i),' iterations'])
                        break
                    end
                end
            end

end
    end 
end
