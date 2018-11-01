function KalmanX2sim(sigmaVis,sigmaSpeed,speedprofile,dt,Ntrial,Fshownoise)
c{1} = 'c';c{2} = 'k';c{3} = 'm';
xmax = 100;
gain = [0.8 1 1.2];%[0.8 1 1.2];%
Tmax = 50;
Rdist_unwrap = zeros(Ntrial,Tmax,numel(gain));
x_unwrap = zeros(Ntrial,Tmax,numel(gain));
Runacc =  zeros(Ntrial,Tmax,numel(gain));
Visacc =  zeros(Ntrial,Tmax,numel(gain));
Runspeed =  speedprofile(1)*ones(Ntrial,Tmax,numel(gain));
Visspeed =  speedprofile(1)*ones(Ntrial,Tmax,numel(gain));
time = zeros(Ntrial,Tmax,numel(gain));
for trial = 1:Ntrial
    for g = 1:numel(gain)
        tt = 2;
        while x_unwrap(trial,tt-1,g) < xmax
            tt = tt+1;
            if tt>100
                keyboard
            end
            time(trial,tt,g) = time(trial,tt-1,g) + dt;
            Runspeed(trial,tt,g) = speedprofile(round(mod(x_unwrap(trial,tt-1,g),xmax))+1);
            Visspeed(trial,tt,g) = gain(g)*speedprofile(round(mod(x_unwrap(trial,tt-1,g),xmax))+1);
            Runacc(trial,tt,g) = (Runspeed(trial,tt,g)-Runspeed(trial,tt-1,g))/dt;
            Visacc(trial,tt,g) = (Visspeed(trial,tt,g)-Visspeed(trial,tt-1,g))/dt;
            Rdist_unwrap(trial,tt,g) = Rdist_unwrap(trial,tt-1,g) + Runspeed(trial,tt,g)*dt;% + Runacc(trial,tt,g)*dt^2/2;%*(time(trial,tt,g)-time(trial,tt-1,g));
            x_unwrap(trial,tt,g) = x_unwrap(trial,tt-1,g) + Visspeed(trial,tt,g)*dt;% + Visacc(trial,tt,g)*dt^2/2;
        end
        Runspeed(trial,1:tt-1,g) = diff(squeeze(Rdist_unwrap(trial,1:tt,g)))/dt;
        Visspeed(trial,1:tt-1,g) = diff(squeeze(x_unwrap(trial,1:tt,g)))/dt;
        Runacc(trial,1:tt-1,g) = diff(squeeze(Runspeed(trial,1:tt,g)))/dt;
        Visacc(trial,1:tt-1,g) = diff(squeeze(Visspeed(trial,1:tt,g)))/dt;
    end
end
x = x_unwrap;
time = time;


duration = (size(x,2)-1) * dt;
Nstate = 2;
figure(1);clf
%% simulate what the Mouse sees over time
for trial = 1:Ntrial
    for g = 1:numel(gain)
        A = [1 dt; 0 1] ; % state transition matrix:  expected trajectory (state prediction)
        B = [dt^2/2; dt]; %input control matrix:  expected effect of the input velocity on the state.
        C = [1 0;0 1];
        
        Q= [0;0]; %initized state--it has two components: [position; velocity] of the Mouse
        Q_estimate = Q;  %x_estimate of initial location estimation of where the Mouse is (what we are updating)
        if nargin <2
            sigmaSpeed = 1; %process noise: the variability in how fast the Mouse is speeding up (stdv of acceleration: meters/sec^2)
            sigmaVis = 2;  %measurement noise: How mask-blinded is the Mouse (stdv of location, in meters)
        end
        Ez = [sigmaVis^2 0; 0 sigmaSpeed^2];% Ez convert the measurement noise (stdv) into covariance matrix
        Ex = sigmaSpeed^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % Ex convert the process noise (stdv) into covariance matrix
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
            Q= A * Q+ B * squeeze(Visacc(trial,round(t/dt)+1,g)) + double(Fshownoise) * MouseSpeed_noise;
%             Q= Q + double(Fshownoise) * MouseSpeed_noise * squeeze(Visacc(trial,round(t/dt)+1,g));
            % Generate what the Mouse sees
            MouseVision_noise = sigmaVis * randn;
            y = C * Q + double(Fshownoise) * MouseVision_noise;
            Q_loc = [Q_loc; Q(1)];
            Q_loc_meas = [Q_loc_meas; y(1)];
            vel = [vel; Q(2)];
            %     acc = [acc; Q(3)];
            %iteratively plot what the Mouse sees
            %     plot(0:dt:t, Q_loc, '-r.')
            %     plot(0:dt:t, Q_loc_meas, '-k.')
            %     axis([0 10 -30 80])
            %     axis tight
            %     hold on
        end
        vel = diff(Q_loc_meas)/dt/gain(g);
        vel = [vel(1);vel];
        
        Q_loc_estimate = []; %  Mouse position estimate
        vel_estimate = []; % Mouse velocity estimate
        Q= [0; 0]; % re-initized state
        P_estimate = P;
        P_mag_estimate = [];
        predic_state = [];
        predic_var = [];
        for t = 1:length(Q_loc)
            % Predict next state of the Mouse with the last state and predicted motion.
            Q_estimate = A * Q_estimate + B * squeeze(Runacc(trial,t,g));
            predic_state = [predic_state; Q_estimate(1)] ;
            %predict next covariance
            P = A * P * A' + Ex;
            predic_var = [predic_var; P] ;
            % predicted Mouse measurement covariance
            % Kalman Gain
            K = P*C'/(C*P*C'+Ez);
            % Update the state estimate.
            Q_estimate = Q_estimate + K * ([Q_loc_meas(t);vel(t)] - C * Q_estimate);
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
        plot(tt,Q_loc,'k',tt,Q_loc_meas/gain(g),['--' c{g}], tt,Q_loc_estimate,['-' c{g} '.']);
        axis tight
        set(gca,'Ylim',[-10 130])
%         subplot(1,numel(gain),2);
%         hold on
%         plot(tt,Q_loc_estimate,['-' c{g} '.'])
    end
end

end