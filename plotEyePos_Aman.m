function plotEyePos_Aman(Foverlaphalves,animalName,varargin)
nAnimal = numel(varargin);
figure('Name',animalName);
for k = 1:nAnimal
    es = varargin{k};    

    eyeXpos = smthInTime(es.eyeXpos, 60, 140, 'same', [], 'median');
    eyeYpos = smthInTime(es.eyeYpos, 60, 140, 'same', [], 'median');
    pupilsize = smthInTime(es.pupilSize, 60, 140, 'same', [], 'median');
    pupilsize = 2*(pupilsize/pi).^0.5;
    pupilDiam = nanmean(pupilsize);
    eyeXpos = eyeXpos/pupilDiam;
    eyeYpos = eyeYpos/pupilDiam;
%     runSpeed = smthInTime(es.smthBallSpd, 60, 140, 'same', [], 'median');
    runSpeed = smthInTime(es.ballspeed, 60, 140, 'same', [], 'median');
    
%     idx = find(es.smthBallSpd>5 & es.contrast>0);
    idx = (es.ballspeed>1 & es.contrast>0);
    
    trialID = unique(es.trialID(idx));
    dx = 1;
    [eyepro,~,~] = fast1Dmap(es.traj(idx),eyeXpos(idx),dx,1,[],false);

    eyeXpostrial = NaN(numel(trialID),size(eyepro,1));
    eyeYpostrial = NaN(numel(trialID),size(eyepro,1));
    pupilsizetrial = NaN(numel(trialID),size(eyepro,1));
    runSpeedtrial = NaN(numel(trialID),size(eyepro,1));
    for tt = 1:numel(trialID)
        [eyepro,~,~] = fast1Dmap(es.traj(idx & es.trialID == trialID(tt)),eyeXpos(idx & es.trialID == trialID(tt)),1,1,[],false);
        eyeXpostrial(tt,1:numel(eyepro)) = eyepro';
        [eyepro,~,~] = fast1Dmap(es.traj(idx & es.trialID == trialID(tt)),eyeYpos(idx & es.trialID == trialID(tt)),1,1,[],false);
        eyeYpostrial(tt,1:numel(eyepro)) = eyepro';
        [eyepro,~,~] = fast1Dmap(es.traj(idx & es.trialID == trialID(tt)),pupilsize(idx & es.trialID == trialID(tt)),1,1,[],false);
        pupilsizetrial(tt,1:numel(eyepro)) = eyepro';
        [spdpro,~,~] = fast1Dmap(es.traj(idx & es.trialID == trialID(tt)),runSpeed(idx & es.trialID == trialID(tt)),1,1,[],false);
        runSpeedtrial(tt,1:numel(eyepro)) = spdpro';
    end
    Xhvalprofile = NaN(1,60);
    Yhvalprofile = NaN(1,60);
    Sizehvalprofile = NaN(1,60);
    Spdhvalprofile = NaN(1,60);
    pval_th = 0.05;
    for ix = 1:60
        eyepos1 = eyeXpostrial(:,ix);
        eyepos2 = eyeXpostrial(:,ix+40);
        [hval,pval] = ttest2(eyepos1(~isnan(eyepos1)),eyepos2(~isnan(eyepos2)),'Alpha',pval_th);
        Xpvalprofile(ix) = pval;
        Xhvalprofile(ix) = hval;

        eyepos1 = eyeYpostrial(:,ix);
        eyepos2 = eyeYpostrial(:,ix+40);
        [hval,pval] = ttest2(eyepos1(~isnan(eyepos1)),eyepos2(~isnan(eyepos2)),'Alpha',pval_th);
        Ypvalprofile(ix) = pval;
        Yhvalprofile(ix) = hval;
        
        eyepos1 = pupilsizetrial(:,ix);
        eyepos2 = pupilsizetrial(:,ix+40);
        [hval,pval] = ttest2(eyepos1(~isnan(eyepos1)),eyepos2(~isnan(eyepos2)),'Alpha',pval_th);
        Sizepvalprofile(ix) = pval;
        Sizehvalprofile(ix) = hval;
        
        spdpos1 = runSpeedtrial(:,ix);
        spdpos2 = runSpeedtrial(:,ix+40);
        [hval,pval] = ttest2(spdpos1(~isnan(spdpos1)),spdpos2(~isnan(spdpos2)),'Alpha',pval_th);
        Spdpvalprofile(ix) = pval;
        Spdhvalprofile(ix) = hval;
    end
    Xhvalprofile(Xhvalprofile == 0) = NaN;
    Yhvalprofile(Yhvalprofile == 0) = NaN;
    Sizehvalprofile(Sizehvalprofile == 0) = NaN;
    Spdhvalprofile(Spdhvalprofile == 0) = NaN;

    [Xeyepro,x,Xeyeprostd] = fast1Dmap(es.traj(idx),eyeXpos(idx),dx,1,[],false);
    [Yeyepro,x,Yeyeprostd] = fast1Dmap(es.traj(idx),eyeYpos(idx),dx,1,[],false);
    [Sizeeyepro,x,Sizeeyeprostd] = fast1Dmap(es.traj(idx),pupilsize(idx),dx,1,[],false);
    [Spdpro,x,Spdprostd] = fast1Dmap(es.traj(idx),runSpeed(idx),dx,1,[],false);
    
    Xeyepro = nanmean(eyeXpostrial,1);
    Yeyepro = nanmean(eyeYpostrial,1);
    Sizeeyepro = nanmean(pupilsizetrial,1);
    Spdpro = nanmean(runSpeedtrial,1);
    
    Xeyeprostd = nanstd(eyeXpostrial,1);
    Yeyeprostd = nanstd(eyeYpostrial,1);
    Sizeeyeprostd = nanstd(pupilsizetrial,1);
    Spdprostd = nanstd(runSpeedtrial,1);
    
    if ~Foverlaphalves
        subplot(4,nAnimal,k);
        ax = gca;
        plot(ax,x,eyeXpostrial(sum(isnan(eyeXpostrial),2) == 0,:),'color',[0.5,0.5,0.5]);
        hold on;
        plot(ax,x,Xeyepro,'k','LineWidth',3);%ciplot(ax,x,Xeyepro,Xeyeprostd,0.1,'k');
        maxdisp = max(abs(eyeXpostrial(:)));%max(max(abs(Xeyepro+Xeyeprostd)),max(abs(Xeyepro-Xeyeprostd)));
        hold on;
        plot(x(10:50),Xhvalprofile(10:50)*(maxdisp*1.2),'k','LineWidth',3);
        plot(x(50:90),Xhvalprofile(10:50)*(maxdisp*1.2),'k','LineWidth',3);
        set(gca,'Xlim',[x(1) x(end)],'Ylim',[-(maxdisp*1.3) (maxdisp*1.3)]);
        if k == 1
            ylabel('Pupil X (A.U.)')
        end
        title([animalName '-' num2str(es.serieslist)])

        subplot(4,nAnimal,nAnimal + k);
        ax = gca;
        plot(ax,x,eyeYpostrial(sum(isnan(eyeYpostrial),2) == 0,:),'color',[0.5,0.5,0.5]);
        hold on;
        plot(ax,x,Yeyepro,'k','LineWidth',3);%ciplot(ax,x,Yeyepro,Yeyeprostd,0.1,'k');
        maxdisp = max(abs(eyeYpostrial(:)));%maxdisp = max(max(abs(Yeyepro+Yeyeprostd)),max(abs(Yeyepro-Yeyeprostd)));
        hold on;
        plot(x(10:50),Yhvalprofile(10:50)*(maxdisp*1.2),'k','LineWidth',3);
        plot(x(50:90),Yhvalprofile(10:50)*(maxdisp*1.2),'k','LineWidth',3);
        set(gca,'Xlim',[x(1) x(end)],'Ylim',[-(maxdisp*1.3) (maxdisp*1.3)]);
        if k == 1
            ylabel('Pupil Y (A.U.)')
        end

        subplot(4,nAnimal,2*nAnimal + k);
        ax = gca;
        plot(ax,x,pupilsizetrial(sum(isnan(pupilsizetrial),2) == 0,:),'color',[0.5,0.5,0.5]);
        hold on;
        plot(ax,x,Sizeeyepro,'k','LineWidth',3);%ciplot(ax,x,Sizeeyepro,Sizeeyeprostd,0.1,'k');
        maxdisp = max((pupilsizetrial(:)));%maxdisp = max(max(abs(Sizeeyepro+Sizeeyeprostd)),max(abs(Sizeeyepro-Sizeeyeprostd)));
        mindisp = min((pupilsizetrial(:)));%mindisp = min(min(abs(Sizeeyepro+Sizeeyeprostd)),min(abs(Sizeeyepro-Sizeeyeprostd)));
        hold on
        plot(x(10:50),Sizehvalprofile(10:50)*(maxdisp*1.01),'k','LineWidth',3);
        plot(x(50:90),Sizehvalprofile(10:50)*(maxdisp*1.01),'k','LineWidth',3);
        set(gca,'Xlim',[x(1) x(end)],'Ylim',[(mindisp*0.95) (maxdisp*1.05)]);
        if k == 1
            ylabel('Pupil Size (A.U.)')
        end

        subplot(4,nAnimal,3*nAnimal + k);
        ax = gca;
        plot(ax,x,runSpeedtrial(sum(isnan(runSpeedtrial),2) == 0,:),'color',[0.5,0.5,0.5]);
        hold on;
        plot(ax,x,Spdpro,'k','LineWidth',3);%ciplot(ax,x,Spdpro,Spdprostd,0.1,'k');
        maxdisp = max(abs(runSpeedtrial(:)));%maxdisp = max(max(abs(Spdpro+Spdprostd)),max(abs(Spdpro-Spdprostd)));
        hold on;
        plot(x(10:50),Spdhvalprofile(10:50)*(maxdisp+0.5),'k','LineWidth',3);
        plot(x(50:90),Spdhvalprofile(10:50)*(maxdisp+0.5),'k','LineWidth',3);
        set(gca,'Xlim',[x(1) x(end)],'Ylim',[0 (maxdisp+1)]);
        if k == 1
            ylabel('Spd (cm/s)')
        end
        xlabel('position (halves, cm)')
    else
        subplot(4,nAnimal,k);
        ax = gca;
        ciplot(ax,x(11:51),Xeyepro(11:51),Xeyeprostd(11:51),0.1,'b');
        hold on;
        ciplot(ax,x(11:51),Xeyepro(51:91),Xeyeprostd(51:91),0.1,'r');
        maxdisp = max(max(abs(Xeyepro+Xeyeprostd)),max(abs(Xeyepro-Xeyeprostd)));
        hold on;
        plot(x(11:51),Xhvalprofile(11:51)*(maxdisp*1.2),'k','LineWidth',3);
        set(gca,'Xlim',[x(11) x(51)],'Ylim',[-(maxdisp*1.3) (maxdisp*1.3)]);
        if k == 1
            ylabel('Pupil X (A.U.)')
        end
        title([animalName '-' num2str(es.serieslist)])

        subplot(4,nAnimal,nAnimal + k);
        ax = gca;
        ciplot(ax,x(11:51),Yeyepro(11:51),Yeyeprostd(11:51),0.1,'b');
        hold on;
        ciplot(ax,x(11:51),Yeyepro(51:91),Yeyeprostd(51:91),0.1,'r');
        maxdisp = max(max(abs(Yeyepro+Yeyeprostd)),max(abs(Yeyepro-Yeyeprostd)));
        hold on;
        plot(x(11:51),Yhvalprofile(11:51)*(maxdisp*1.2),'k','LineWidth',3);
        set(gca,'Xlim',[x(11) x(51)],'Ylim',[-(maxdisp*1.3) (maxdisp*1.3)]);
        if k == 1
            ylabel('Pupil Y (A.U.)')
        end

        subplot(4,nAnimal,2*nAnimal + k);
        ax = gca;
        ciplot(ax,x(11:51),Sizeeyepro(11:51),Sizeeyeprostd(11:51),0.1,'b');
        hold on;
        ciplot(ax,x(11:51),Sizeeyepro(51:91),Sizeeyeprostd(51:91),0.1,'r');
        maxdisp = max(max(abs(Sizeeyepro+Sizeeyeprostd)),max(abs(Sizeeyepro-Sizeeyeprostd)));
        mindisp = min(min(abs(Sizeeyepro+Sizeeyeprostd)),min(abs(Sizeeyepro-Sizeeyeprostd)));
        hold on;
        plot(x(11:51),Sizehvalprofile(11:51)*(maxdisp*1.01),'k','LineWidth',3);
        set(gca,'Xlim',[x(11) x(51)],'Ylim',[(mindisp*0.95) (maxdisp*1.05)]);
        if k == 1
            ylabel('Pupil Size (A.U.)')
        end

        subplot(4,nAnimal,3*nAnimal + k);
        ax = gca;
        ciplot(ax,x(11:51),Spdpro(11:51),Spdprostd(11:51),0.1,'b');
        hold on;
        ciplot(ax,x(11:51),Spdpro(51:91),Spdprostd(51:91),0.1,'r');
        maxdisp = max(max(abs(Spdpro+Spdprostd)),max(abs(Spdpro-Spdprostd)));
        hold on;
        plot(x(11:51),Spdhvalprofile(11:51)*(maxdisp*1.01),'k','LineWidth',3);
        set(gca,'Xlim',[x(11) x(51)],'Ylim',[0 (maxdisp*1.05)]);
        if k == 1
            ylabel('Spd (cm/s)')
        end
        xlabel('position (halves, cm)')
    end
end