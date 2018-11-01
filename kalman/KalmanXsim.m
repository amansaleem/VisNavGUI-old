function KalmanXsim(sigmaVis,sigmaSpeed,speedprofile,dt,Ntrial,Fshownoise)
c{1} = 'c';c{2} = 'k';c{3} = 'm';
xmax = 100;
gain = [0.8 1 1.2];%[0.8 1 1.2];%
Tmax = 50;
Rdist_unwrap = zeros(Ntrial,Tmax,numel(gain));
x_unwrap = zeros(Ntrial,Tmax,numel(gain));
Runspeed =  zeros(Ntrial,Tmax,numel(gain));
Visspeed =  zeros(Ntrial,Tmax,numel(gain));
time = zeros(Ntrial,Tmax,numel(gain));
for trial = 1:Ntrial
    for g = 1:numel(gain)
        tt = 2;
        while x_unwrap(trial,tt-1,g) < xmax
            tt = tt+1;
            time(trial,tt,g) = time(trial,tt-1,g) + dt;
            Runspeed(trial,tt,g) = speedprofile(round(mod(x_unwrap(trial,tt-1,g),xmax))+1);
            Visspeed(trial,tt,g) = gain(g)*speedprofile(round(mod(x_unwrap(trial,tt-1,g),xmax))+1);
            
            Rdist_unwrap(trial,tt,g) = Rdist_unwrap(trial,tt-1,g) + Runspeed(trial,tt,g)*dt;%*(time(trial,tt,g)-time(trial,tt-1,g));
            x_unwrap(trial,tt,g) = x_unwrap(trial,tt-1,g) + Visspeed(trial,tt,g)*dt;
        end
    end
end
x = x_unwrap;
time = time;

duration = (size(x,2)-1) * dt;
Nstate = 1;
figure(1);clf
Q_loc_estimate_all = cell(1,3);
Q_loc_all = cell(1,3);
%% simulate what the Mouse sees over time
for trial = 1:Ntrial
    for g = 1:numel(gain)
        
        A = [1] ; % state transition matrix:  expected trajectory (state prediction)
        B = [dt]; %input control matrix:  expected effect of the input velocity on the state.
        C = [1];
        
        Q= [0]; %initized state--it has two components: [position; velocity] of the Mouse
        Q_estimate = Q;  %x_estimate of initial location estimation of where the Mouse is (what we are updating)
        if nargin <2
            sigmaSpeed = 1; %process noise: the variability in how fast the Mouse is speeding up (stdv of acceleration: meters/sec^2)
            sigmaVis = 2;  %measurement noise: How mask-blinded is the Mouse (stdv of location, in meters)
        end
        Ez = (sigmaVis)^2;% Ez convert the measurement noise (stdv) into covariance matrix
        Ex = sigmaSpeed^2 * [dt^2]; % Ex convert the process noise (stdv) into covariance matrix
        P = Ex; % estimate of initial Mouse position variance (covariance matrix)
        
        %% initize result variables
        % Initialize for speed
        Q_loc = []; % ACTUAL Mouse flight path
        vel = []; % ACTUAL Mouse velocity
        acc = [];
        Q_loc_meas = []; % Mouse path that the Mouse sees


        
        duration = (find(squeeze(x(trial,:,g)) > 0,1,'last')-1)*dt;
        for t = 0 : dt: duration
            
            % Generate the Mouse flight
            MouseSpeed_noise = sigmaSpeed * [(dt^2)*randn];
            Q= A * Q+ B * squeeze(Visspeed(trial,round(t/dt)+1,g)) + double(Fshownoise) * MouseSpeed_noise;
            % Generate what the Mouse sees
            MouseVision_noise = sigmaVis * randn;
            y = C * Q + double(Fshownoise) * MouseVision_noise;
            Q_loc = [Q_loc; Q(1)];
            Q_loc_meas = [Q_loc_meas; y];
            %     vel = [vel; Q(2)];
            %     acc = [acc; Q(3)];
            %iteratively plot what the Mouse sees
            %     plot(0:dt:t, Q_loc, '-r.')
            %     plot(0:dt:t, Q_loc_meas, '-k.')
            %     axis([0 10 -30 80])
            %     axis tight
            %     hold on
        end
        
        Q_loc_estimate = []; %  Mouse position estimate
        vel_estimate = []; % Mouse velocity estimate
        Q= [0; 0]; % re-initized state
        P_estimate = P;
        P_mag_estimate = [];
        predic_state = [];
        predic_var = [];
        for t = 1:length(Q_loc)
            % Predict next state of the Mouse with the last state and predicted motion.
            Q_estimate = A * Q_estimate + B * squeeze(Runspeed(trial,t,g));
            predic_state = [predic_state; Q_estimate(1)] ;
            %predict next covariance
            P = A * P * A' + Ex;
            predic_var = [predic_var; P] ;
            % predicted Mouse measurement covariance
            % Kalman Gain
            if (Q_loc_meas(t)<10) || (Q_loc_meas(t)>35 && Q_loc_meas(t)<45) || (Q_loc_meas(t)>55 && Q_loc_meas(t)<65) || (Q_loc_meas(t)>90)
                Ezscaling = 500;
            else
                Ezscaling = 1;
            end
            K = P*C'/(C*P*C'+Ez/Ezscaling);
            % Update the state estimate.
            Q_estimate = Q_estimate + K * (Q_loc_meas(t) - C * Q_estimate);
%             if t > 1 && t < length(Q_loc)
%                 Runspeed(trial,t+1,g) = 0.5*(Runspeed(trial,t+1,g) + (Q_estimate(1) - Q_loc_estimate(t-1))/dt);
%             end
            % update covariance estimation.
            P =  (eye(Nstate)-K*C)*P;
            %Store for plotting
            Q_loc_estimate = [Q_loc_estimate; Q_estimate(1)];
            %     vel_estimate = [vel_estimate; Q_estimate(2)];
            P_mag_estimate = [P_mag_estimate; P(1)];
        end
        % figure(1)
        % hold on
        % plot(0:dt:duration, vel_estimate, '-g.');
        % Plot the results
        subplot(1,numel(gain),g);
        tt = 0 : dt : duration;
        plot(tt,Q_loc,'k',tt,Q_loc_meas/gain(g),['--' c{g}], tt,Q_loc_estimate,['-' c{g} '.'],tt,Q_loc_estimate-Q_loc,c{g},tt,tt*0,'k');
        axis tight
        set(gca,'Ylim',[-10 130])
        subplot(1,numel(gain),2);
        hold on
        plot(tt,Q_loc_estimate,['-' c{g} '.'])
        Q_loc_all{g} = [Q_loc_all{g};Q_loc];
        Q_loc_estimate_all{g} = [Q_loc_estimate_all{g};Q_loc_estimate];
    end
end
figure;
lambdaSmooth = 2;
for g = 1:3
    [F, ~, ~, ~] = smoothhist2D_corrected([Q_loc_all{g} Q_loc_estimate_all{g}], lambdaSmooth, [100 100], 1:100, 1:100, true, true);
    subplot(2,3,g);imagesc(F)
    Xpred{g} = getCircularAverage(F,1,0.1);
    set(gca,'Ydir','normal','Clim',[0 4])
    hold on;plot([1 100], [1 100],'k');
    subplot(2,3,4);
    hold on
    plot((Xpred{g})'-(1:100),c{g});
    set(gca,'Ylim',[-20 20])
end
for g = [1 3]
    [Xpred_reg,x] = Xpredcorrection(Xpred{g},Xpred{2});
    subplot(2,3,6);
    hold on
    plot(x,Xpred_reg,c{g});
    set(gca,'Ylim',[-20 20])
end
end

function [predave_out,x] = Xpredcorrection(predave,predave_ref)
predave_ref(predave_ref - (1:100)' > 50) = predave_ref(predave_ref - (1:100)' > 50) - 100;
predave_ref(predave_ref - (1:100)' < -50) = predave_ref(predave_ref - (1:100)' < -50) + 100;
predave(predave - (1:100)' > 50) = predave(predave - (1:100)' > 50) - 100;
predave(predave - (1:100)' < -50) = predave(predave - (1:100)' < -50) + 100;

predave_interp = interp1(linspace(0,numel(predave),numel(predave)+1), [predave(1);predave], 0:0.1:99.9);predave_interp(isnan(predave_interp)) = 0;
predaveref_interp = interp1(linspace(0,numel(predave_ref),numel(predave_ref)+1), [predave_ref(1);predave_ref], 0:0.1:99.9);predaveref_interp(isnan(predaveref_interp)) = 0;
predave_out = zeros(size(predave_interp));
x = 0:0.1:99.9;
predaveref_interp = [predaveref_interp(1:end-1)-100 predaveref_interp predaveref_interp(2:end)+100];
xrep = [x(1:end-1)-100 x x(2:end)+100];
for i = 1:numel(x)
    try
    idxmatch = min(find(abs(predaveref_interp-predave_interp(i)) <= min(abs(predaveref_interp-predave_interp(i)))));
    idxmatch = idxmatch(abs(idxmatch-(numel(x)-1) - i) == min(abs(idxmatch-(numel(x)-1) - i)));
    predave_out(i) = xrep(round(idxmatch)) - x(i);
    catch
        keyboard
    end
end
% predave_out(predave_out>50) = predave_out(predave_out>50) - 100;
% predave_out(predave_out<-50) = predave_out(predave_out<-50) + 100;
end