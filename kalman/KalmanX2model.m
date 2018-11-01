function Xmodel = KalmanX2model(Xactual,sigmaVis,sigmaSpeed,sigmaAcc,latcorrection,es,Tsmth_win, wintype,numBins)
sampleRate = mean(1./es.sampleSize);
dt = mean(es.sampleSize);

Xactual = circshift(Xactual,[round(latcorrection/(1000/sampleRate)) 0]);
if sum(sigmaVis) == 0 && sigmaSpeed == 0
    Xmodel = Xactual;
else
%     Xtemp = unwrap(Xactual/max(Xactual)*2*pi)*max(Xactual)/(2*pi);
%     Xtemp = smthInTime(Xtemp, sampleRate, Tsmth_win, 'same', [], wintype);
%     Xtemp = mod(Xtemp,max(Xactual));
    Xtemp = Xactual;

    trialID = unique(es.trialID);
    Xmodel = zeros(size(Xtemp));
    maxX = max(Xtemp);
    for trial = 1:numel(trialID)
        trialidx = min(find(es.trialID == trialID(trial)),numel(Xtemp));
        Xtrial = Xtemp(trialidx);
        Xtrial = unwrap(Xtrial/maxX*2*pi)*maxX/(2*pi);
        gaintrial = es.gain(trialidx)/mode(es.gain);

        Runspeed = [0;diff(Xtrial)/dt]./(gaintrial);
        Runacc = [0;diff(Runspeed)/dt];
%         Runacc = smooth(Runacc,10);
%         Runacc = circshift(Runacc,0);
        Nstate = 2;
        A = [1 dt; 0 1] ; % state transition matrix:  expected trajectory (state prediction)
        B = [dt^2/2; dt]; %input control matrix:  expected effect of the input velocity on the state.
        C = [1 0;0 1];

        X= [Xtrial(1); Runspeed(1)]; %initized state--it has two components: [position; velocity] of the Mouse
        X_estimate = X;  %x_estimate of initial location estimation of where the Mouse is (what we are updating)
        if nargin <2
            sigmaSpeed = 1; %process noise: the variability in how fast the Mouse is speeding up (stdv of acceleration: meters/sec^2)
            sigmaVis = 2;  %measurement noise: How mask-blinded is the Mouse (stdv of location, in meters)
        end
        nSigmaVis = numel(sigmaVis);
        if nSigmaVis == 1
            Ez = [sigmaVis^2 0; 0 sigmaSpeed^2];% Ez convert the measurement noise (stdv) into covariance matrix
        end
        Ex = sigmaAcc^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % Ex convert the process noise (stdv) into covariance matrix
        P = Ex; % estimate of initial Mouse position variance (covariance matrix)

        Xmodeltrial = []; %  Mouse position estimate
        vel_estimate = []; % Mouse velocity estimate
        P_estimate = P;
        P_mag_estimate = [];
        predic_state = [];
        predic_var = [];
        for t = 1:length(Xtrial)
            % Predict next state of the Mouse with the last state and predicted motion.
            X_estimate = A * X_estimate + B * Runacc(t);
            predic_state = [predic_state; X_estimate(1)] ;
            %predict next covariance
            P = A * P * A' + Ex;
            predic_var = [predic_var; P] ;
            % predicted Mouse measurement covariance
            % Kalman Gain
%             if Xtrial(t)>0 && Xtrial(t)<35
%                 Ez(1,1) = 2;
%             else
%                 Ez(1,1) = 0.5;
%             end
            Ez =  [sigmaVis(min(trialidx(t),nSigmaVis))^2 0; 0 sigmaSpeed^2];
            K = P*C'/(C*P*C'+Ez);
            % Update the state estimate.
    %         if mod(t,2) == 0
            X_estimate = X_estimate + K * ([Xtrial(t);Runspeed(t)] - C * X_estimate);        
    %         end
            % update covariance estimation.
            P =  (eye(Nstate)-K*C)*P;
            %Store for plotting
            Xmodeltrial = [Xmodeltrial; X_estimate(1)];
            %     vel_estimate = [vel_estimate; X_estimate(2)];
            P_mag_estimate = [P_mag_estimate; P(1)];
        end
        Xmodel(trialidx) =  mod(Xmodeltrial,maxX);
    end
end
Xmodel = unwrap(Xmodel/max(Xactual)*2*pi)*max(Xactual)/(2*pi);
Xmodel = smthInTime(Xmodel, sampleRate, Tsmth_win, 'same', [], wintype);
Xmodel = mod(Xmodel,max(Xactual));

if numBins > 0
    [Xmodel, ~] = normalise1var(Xmodel, numBins);
end

end