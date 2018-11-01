function Xmodel = KalmanXmodel(Xactual,sigmaVis,sigmaSpeed,latcorrection,es,Tsmth_win, wintype,numBins)
sampleRate = mean(1./es.sampleSize);
dt = mean(es.sampleSize);

Xactual = circshift(Xactual,[round(latcorrection/(1000/sampleRate)) 0]);
if sigmaVis == 0 && sigmaSpeed == 0
    Xmodel = Xactual;
else
    
%     Xtemp = unwrap(Xactual/max(Xactual)*2*pi)*max(Xactual)/(2*pi);
%     Xtemp = smthInTime(Xtemp, sampleRate, Tsmth_win, 'same', [], wintype);
%     Xtemp = mod(Xtemp,max(Xactual));
%     Xtemp = Xactual;
    
    trialID = unique(es.trialID);
    Xmodel = zeros(size(Xactual));
    maxX = max(Xactual);
    for trial = 1:numel(trialID)
        trialidx = min(find(es.trialID == trialID(trial)),numel(Xactual));%1:numel(Xtemp);%
        Xtrial = Xactual(trialidx);
        Xtrial = unwrap(Xtrial/maxX*2*pi)*maxX/(2*pi);
        gaintrial = es.gain(trialidx)/mode(es.gain);
        
        Runspeed = [0;diff(Xtrial)/dt]./(gaintrial);
        Runspeed(1) = Runspeed(2);
%         Runspeed = circshift(Runspeed,5);
        Nstate = 1;
        A = [1] ; % state transition matrix:  expected trajectory (state prediction)
        B = [dt]; %input control matrix:  expected effect of the input velocity on the state.
        C = [1];
        
        X= [Xtrial(1)]; %initized state--it has two components: [position; velocity] of the Mouse
        X_estimate = X;  %x_estimate of initial location estimation of where the Mouse is (what we are updating)
        if nargin <2
            sigmaSpeed = 1; %process noise: the variability in how fast the Mouse is speeding up (stdv of acceleration: meters/sec^2)
            sigmaVis = 2;  %measurement noise: How mask-blinded is the Mouse (stdv of location, in meters)
        end
        Ez = (sigmaVis)^2*ones(size(Xtrial));% Ez convert the measurement noise (stdv) into covariance matrix
        Ex = sigmaSpeed^2 * [dt^2]; % Ex convert the process noise (stdv) into covariance matrix
        P = Ex; % estimate of initial Mouse position variance (covariance matrix)
        
        Xmodeltrial = []; %  Mouse position estimate
        vel_estimate = []; % Mouse velocity estimate
        Q= [0; 0]; % re-initized state
        P_estimate = P;
        P_mag_estimate = [];
        predic_state = [];
        predic_var = [];
        K = 1;
        for t = 1:length(Xtrial)
            % Predict next state of the Mouse with the last state and predicted motion.
            X_estimate = A * X_estimate + B * Runspeed(t);
%             if es.reward(t) == 2
%                 X_estimate = 95;
%             end
            predic_state = [predic_state; X_estimate(1)] ;
            %predict next covariance
            P = A * P * A' + Ex;
            predic_var = [predic_var; P] ;
            % predicted Mouse measurement covariance
            % Kalman Gain
%             if es.reward(t) == 3
%                  K = P*C'/(C*P*C'+0);
%             else
                K = P*C'/(C*P*C'+Ez(t));
%             end
            % Update the state estimate.
            X_estimate = X_estimate + K * (Xtrial(t) - C * X_estimate);
            
            % update covariance estimation.
            P =  (eye(Nstate)-K*C)*P;
            %Store for plotting
            Xmodeltrial = [Xmodeltrial; X_estimate(1)];
            %     vel_estimate = [vel_estimate; X_estimate(2)];
            P_mag_estimate = [P_mag_estimate; P(1)];
        end
%         Xmodel = Xmodeltrial;
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