function [Xfit,sigmaVis,sigmaSpeed] = KalmanXfit(Xactual,Xdec,latcorrection,es,Tsmth_win, wintype,numBins)
sampleRate = mean(1./es.sampleSize);
dt = mean(es.sampleSize);

Xactual = circshift(Xactual,[round(latcorrection/(1000/sampleRate)) 0]);
sigmaVis0 = 1;
sigmaSpeed0 = 1;
    
%     Xtemp = unwrap(Xactual/max(Xactual)*2*pi)*max(Xactual)/(2*pi);
%     Xtemp = smthInTime(Xtemp, sampleRate, Tsmth_win, 'same', [], wintype);
%     Xtemp = mod(Xtemp,max(Xactual));
    
trialID = unique(es.trialID);
Xfit = zeros(size(Xactual));
sigmaVis = zeros(size(Xactual));
sigmaSpeed = zeros(size(Xactual));
maxX = max(Xactual);
for trial = 1:numel(trialID)
    trialidx = min(find(es.trialID == trialID(trial)),numel(Xactual));%1:numel(Xtemp);%
    Xtrial = Xactual(trialidx);
%     Xtrial = unwrap(Xtrial/maxX*2*pi)*maxX/(2*pi);
    Xdectrial = Xdec(trialidx);
%     Xdectrial = unwrap(Xdectrial/maxX*2*pi)*maxX/(2*pi);
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
    Ez = (sigmaVis0)^2;% Ez convert the measurement noise (stdv) into covariance matrix
    Ex = sigmaSpeed0^2 * [dt^2]; % Ex convert the process noise (stdv) into covariance matrix
    P = Ex; % estimate of initial Mouse position variance (covariance matrix)
    
    Xfittrial = []; %  Mouse position estimate
    sigmaVistrial = [];
    sigmaSpeedtrial = [];
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
        K = P*C'/(C*P*C'+Ez);
        %             end
        % Update the state estimate.
        if (Xtrial(t) - C * X_estimate) ~= 0
            K = (Xdectrial(t) - X_estimate)/(Xtrial(t) - C * X_estimate);
            Ez = 1/(K/P*C') - (C*P*C');
        else
            X_estimate = X_estimate + K * (Xtrial(t) - C * X_estimate);
        end
        try
        sigmaVistrial = [sigmaVistrial; Ez];
        catch
            keyboard
        end
        sigmaSpeedtrial = [sigmaSpeedtrial ; 1];
        % update covariance estimation.
        P =  (eye(Nstate)-K*C)*P;
        %Store for plotting
        Xfittrial = [Xfittrial; X_estimate(1)];
        %     vel_estimate = [vel_estimate; X_estimate(2)];
        P_mag_estimate = [P_mag_estimate; P(1)];
    end
    %         Xfit = Xmodeltrial;
    Xfit(trialidx) =  mod(Xfittrial,maxX);
    sigmaVis(trialidx) =  sigmaVistrial;
    sigmaSpeed(trialidx) =  sigmaSpeedtrial;
end

Xfit = unwrap(Xfit/max(Xactual)*2*pi)*max(Xactual)/(2*pi);
Xfit = smthInTime(Xfit, sampleRate, Tsmth_win, 'same', [], wintype);
Xfit = mod(Xfit,max(Xactual));

if numBins > 0
    [Xfit, ~] = normalise1var(Xfit, numBins);
end

end