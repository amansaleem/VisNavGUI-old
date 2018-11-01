function reslat = BatchLatEstim
SetDirs;
expt = getExperimentList;
nanimal = numel(expt);
maxtol = 1;
ampth = 0;
reslat.Latopt = cell(2,3);
reslat.MeanPostXMaxLick = cell(2,3);
reslat.MeanPostXAveLick = cell(2,3);
reslat.LickXmax = cell(2,3);
reslat.LickXave = cell(2,3);
datadir = 'D:\DATA\batch';
filesuffix = '_Twin250_Xwin4_spdth5_3speedbins_3eyebins_0thetabins_goodonly_decoderParams_latency';
figure;
for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        filename = [datadir filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed' filesep animalname '_' num2str(expt(ianimal).series{iseries}) filesuffix '.mat'];
        if exist(filename,'file')
            for iprobe = 1:2
                if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                    S = load(filename);
                    for g = 1:3
                        if ((iprobe == 1 && expt(ianimal).goodCA1dec{g}{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1dec{g}{iseries} == 1))
                            imax{g} = getCircularAverage(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,g,:))',ampth,maxtol);
                            %                     imax{g} = getCircularAverage(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,1,:))',ampth,maxtol);
                            %                     imax{g} = getCircularAverage(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,3,:))',ampth,maxtol);
                        else
                            imax{g} = NaN;
                        end
                    end
%                     [~,imax_med] = max(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,2,:)),[],2);
%                     [~,imax_low] = max(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,1,:)),[],2);
%                     [~,imax_high] = max(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,3,:)),[],2);
                    shift_low = imax{1} - imax{2};
                    shift_high = imax{3} - imax{2};
                    reslat.Latopt{iprobe,1} = [reslat.Latopt{iprobe,1}  mean(S.latparams.latcorrectionlist(find(abs(shift_low) == min(abs(shift_low)),1,'last')))];
                    reslat.Latopt{iprobe,3} = [reslat.Latopt{iprobe,3}  mean(S.latparams.latcorrectionlist(find(abs(shift_high) == min(abs(shift_high)),1,'last')))];
                    reslat.Latopt{iprobe,2} = [reslat.Latopt{iprobe,2}  S.latEstim.OptLatency];
%                     for g = 1:3
%                         imax = getCircularAverage(squeeze(S.latEstim.ParamsLickdistri(iprobe,S.latparams.latcorrectionlist == 0,1,1,1,1,g,:)),0,maxtol);
%                         reslat.LickXave{iprobe,g} = [reslat.LickXave{iprobe,g} S.latparams.xParams(floor(imax))];
%                         imax = getCircularAverage(squeeze(S.latEstim.ParamsLickdistri(iprobe,S.latparams.latcorrectionlist == 0,1,1,1,1,g,:)),0,0.1);
%                         reslat.LickXmax{iprobe,g} = [reslat.LickXmax{iprobe,g} S.latparams.xParams(floor(imax))];
%                         imax = getCircularAverage(squeeze(S.latEstim.ParamsMeanPostbeforelick(iprobe,S.latparams.latcorrectionlist == 0,1,1,1,1,g,:,:))',0,maxtol);
%                         reslat.MeanPostXAveLick{iprobe,g} = [reslat.MeanPostXAveLick{iprobe,g} S.latparams.xParams(floor(imax))'];
%                         imax = getCircularAverage(squeeze(S.latEstim.ParamsMeanPostbeforelick(iprobe,S.latparams.latcorrectionlist == 0,1,1,1,1,g,:,:))',0,0.1);
%                         reslat.MeanPostXMaxLick{iprobe,g} = [reslat.MeanPostXMaxLick{iprobe,g} S.latparams.xParams(floor(imax))'];
%                     end
                else
                    for g = 1:3
                        reslat.Latopt{iprobe,g} = [reslat.Latopt{iprobe,g}  NaN];
                        reslat.LickXmax{iprobe,g} = [reslat.LickXmax{iprobe,g} NaN];
                        reslat.LickXave{iprobe,g} = [reslat.LickXave{iprobe,g} NaN];
                        reslat.MeanPostXMaxLick{iprobe,g} = [reslat.MeanPostXMaxLick{iprobe,g} NaN(60,1)];
                        reslat.MeanPostXAveLick{iprobe,g} = [reslat.MeanPostXAveLick{iprobe,g} NaN(60,1)];
                    end
                end
            end
        end
    end
end
end