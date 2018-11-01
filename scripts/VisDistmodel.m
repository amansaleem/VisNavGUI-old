function [x_hat, x] = VisDistmodel(alpha,delay,speed,gain,dt,Ntrial)
    figure;
    c{1} = '+c';c{2} = '+b';c{3} = '+m';
    xmax = 100;
    gain = [0.8 1 1.2];%
    Tmax = 100;
    Rdist_unwrap = zeros(Ntrial,Tmax,numel(gain));
    x_unwrap = zeros(Ntrial,Tmax,numel(gain));
    time = zeros(Ntrial,Tmax,numel(gain));
    for trial = 1:Ntrial 
        for g = 1:numel(gain)
            tt = 2;
            while x_unwrap(trial,tt-1,g) < xmax
                tt = tt+1;
                time(trial,tt,g) = time(trial,tt-1,g) + dt;
                Rdist_unwrap(trial,tt,g) = Rdist_unwrap(trial,tt-1,g) + speed(round(mod(x_unwrap(trial,tt-1,g),xmax))+1)*dt;%*(time(trial,tt,g)-time(trial,tt-1,g));
                x_unwrap(trial,tt,g) = Rdist_unwrap(trial,tt,g)*gain(g);                
            end
        end
    end
    time = time;
%     xmax = round(max(x_unwrap(:)));
    
%     Rdist_unwrap = [0;cumsum(speed(2:end).*diff(time))];%cumsum(speed.*[diff(time) time(2)-time(1)]);
%     Rdist_unwrap = Rdist_unwrap(:);
%     x_unwrap = Rdist_unwrap(:)*gain;
    Rdist = Rdist_unwrap;%*ones(size(gain));%zeros(size(x_unwrap));
    x = x_unwrap;%zeros(size(x_unwrap));
%     for g = 1:numel(gain)
%         Rdist(:,g) = mod(Rdist_unwrap,xmax./gain(g)); 
%         x(:,g) = Rdist(:,g)*gain(g);
%     end
    
%     timediff = dt;
%     dt = round(mean(timediff>0));
    Ntau = size(x,2);
    delay = 1+round(delay/dt);%1/dt;%+round(delay/dt);
    Xmemory = zeros(1,delay+Ntau);
    Xmemory(delay+1:end) = (alpha).^((0:(Ntau-1))); 
    Xmemory = Xmemory(:);
        
    subplot(2,3,1); 
    for g = 1:numel(gain)
        subplot(2,3,g);
        hold(gca,'on');
        plot(time(:,:,g),x(:,:,g),'.k');
    end
    
    x_hat = cell(3,1);
    for trial = 1:Ntrial
        for g = 1:numel(gain)
            disttemp = squeeze(Rdist_unwrap(trial,:,g));%squeeze(Rdist_unwrap(trial,:,g));
            triallength = find(disttemp > 0,1,'last');
            disttemp = disttemp(1:triallength);
            disttemp = disttemp(:);
            dist = [diff(disttemp);disttemp(end)-disttemp(end-1)];
            dist_model = conv(Xmemory,dist,'full');
            dist_model = dist_model(1:numel(dist));
            
            xtemp = squeeze(x_unwrap(trial,:,g));
            xtemp = xtemp(1:triallength);
            xtemp = xtemp(:);
            x_model = conv(Xmemory,xtemp,'full');
            x_model = x_model(1:numel(xtemp))*(1-alpha);%*dt;
            
            x_hat{trial,g} = x_model + dist_model*(alpha>0);
            
            subplot(2,3,g);
            hold(gca,'on');
            %         plot(x(idx_noedge,g),dist(idx_noedge),'k')
            plot(time(trial,1:triallength,g),x_hat{trial,g},c{g}(2))
            subplot(2,3,4);
            hold(gca,'on');
            plot(xtemp,x_hat{trial,g} - xtemp,c{g}(2))
        end
    end
    oversampling = 10;
    x_hatmean = zeros(xmax*oversampling,numel(gain));
    
    for g = 1:numel(gain)
        nXX = zeros(1,xmax*oversampling);
        for xx = 2:xmax*oversampling
            for trial = 1:Ntrial
                x_hatmean(xx,g) = x_hatmean(xx,g) + sum(x_hat{trial,g}(squeeze(round(x_unwrap(trial,:,g)*oversampling)) == xx));
                try
                nXX(xx) = nXX(xx) + sum(squeeze(round(x_unwrap(trial,:,g)*oversampling)) == xx);
                catch
                    keyboard
                end
            end
            x_hatmean(xx,g) = x_hatmean(xx,g)/nXX(xx);
        end 
        goodval = find(~isnan(x_hatmean(:,g)));
        x_hatmean(:,g) = interp1(goodval/oversampling,x_hatmean(goodval,g),linspace(0,99,100*oversampling),'linear');
        
        subplot(2,3,5);
        hold(gca,'on'); 
        plot(x_hatmean(:,g),c{g})
    end
    
    for trial = 1:Ntrial
        predave1 = x_hatmean(:,1);
        predave2 = x_hatmean(:,2);
        predave3 = x_hatmean(:,3);
        
        predave1_interp = interp1(linspace(0,100,numel(predave1)+1), [predave1(1);predave1], 0:0.1:99.9);predave1_interp(isnan(predave1_interp)) = 0;
        predave2_interp = interp1(linspace(0,100,numel(predave2)+1), [predave2(1);predave2], 0:0.1:99.9);predave2_interp(isnan(predave2_interp)) = 0;
        predave3_interp = interp1(linspace(0,100,numel(predave3)+1), [predave3(1);predave3], 0:0.1:99.9);predave3_interp(isnan(predave3_interp)) = 0;
        predave_1 = zeros(size(predave1_interp));
        predave_2 = zeros(size(predave2_interp));
        predave_3 = zeros(size(predave3_interp));
        x = 0:0.1:99.9;
        for i = 1:numel(x)
            predave_1(i) = x(round(mean(find(abs(predave2_interp-predave1_interp(i)) <= min(abs(predave2_interp-predave1_interp(i))))))) - x(i);
            predave_2(i) = x(round(mean(find(abs(predave2_interp-predave2_interp(i)) <= min(abs(predave2_interp-predave2_interp(i))))))) - x(i);
            predave_3(i) = x(round(mean(find(abs(predave2_interp-predave3_interp(i)) <= min(abs(predave2_interp-predave3_interp(i))))))) - x(i);
        end
        
        subplot(2,3,6);
        hold(gca,'on'); 
        plot(x,predave_1,c{1}(2));
        plot(x,predave_2,c{2}(2))
        plot(x,predave_3,c{3}(2))
    end
    
    
%     for g = 1:numel(gain)
%         subplot(2,3,6);
%         hold(gca,'on');         
%         plot(x_hatmean(:,g)-x_hatmean(:,2),c{g})
% %         plot(x(x(:,g)<xmax,g),x_hat{g}(x(:,g)<xmax),c{g}(2));
% %         plot(mod(x(:,g),xmax),mod(x(:,g),xmax) + (x_hat{g}(:)-x(:,g)),c{g});
% %         plot(x(x(:,g)<xmax,g),x_hat{g}(x(:,g)<xmax)-x(x(:,g)<xmax,g) - (x_hat{2}(x(:,g)<xmax)-x(x(:,g)<xmax,2)),c{g});
%     end        
end