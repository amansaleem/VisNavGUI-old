function es_spk = process_BusseData
% % •	V1, c, gray: BL6_0188,  s7, e7
% % •	V1, a, gray: BL6_0188,  s7, e10
% % % •	V1, c, black: 'Ntsr1-Cre_0039',  s9, e3
% % % •	V1, c, black: BL6_0188,  s7, e6
% % % •	V1, c, black: 'Ntsr1-Cre_0039',  s9, e4
% % % •	V1, c, black: 'Ntsr1-Cre_0039',  s10, e2
% % % •	V1, a, black: 'Ntsr1-Cre_0039',  s9, e6
% % % •	V1, a, black: BL6_0188,  s7, e13
% % % •	V1, a, black: 'Ntsr1-Cre_0039',  s9, e7
% % % •	V1, a, black: 'Ntsr1-Cre_0039',  s10, e3
% % •	LGN, c, gray: BL6_0191,  s5, e2
% % •	LGN, c, gray: BL6_0191,  s6, e4
% % •	LGN, a, gray: BL6_0191,  s5, e6
% % •	LGN, a, gray: BL6_0191,  s6, e9
% % •	LGN, c, black: BL6_0191,  s5, e4
% % •	LGN, c, black: BL6_0191,  s6, e5
% % •	LGN, a, black: BL6_0191,  s5, e7
% % •	LGN, a, black: BL6_0191,  s6, e10
global d
dirName  = 'C:\Users\aman\Dropbox\Work\Data\V1 gamma\Data';

if isempty(d)
    load([dirName filesep 'forAman']);     
end
% % % LGN Black, dilated
% % doOneFile(dirName, d.lgn_a_black(1))
% % set(gcf,'name','LGN, Atropine Black 1')
% % 
% % doOneFile(dirName, d.lgn_a_black(2))
% % set(gcf,'name','LGN, Atropine Black 2')
% % 
% % % LGN Black, control
% % doOneFile(dirName, d.lgn_c_black(1))
% % set(gcf,'name','LGN, Control Black 1')
% % 
% % doOneFile(dirName, d.lgn_c_black(2))
% % set(gcf,'name','LGN, Control Black 2')

% % LGN Gray, dilated
% doOneFile(dirName, d.lgn_a_gray(1), 20, [10 20 22 25])
% set(gcf,'name','LGN, Atropine Gray 1')
% 
% es_spk.a = doOneFile(dirName, d.lgn_a_gray(2), 20, [7 9])
% set(gcf,'name','LGN, Atropine Gray 2')

% LGN Gray, control
es_spk.c = doOneFile(dirName, d.lgn_c_gray(2), 20);%[7 9]);
set(gcf,'name','LGN, control Gray 2')

doOneFile(dirName, d.lgn_c_gray(1), 20);%[10 20 22 25])
set(gcf,'name','LGN, control Gray 1')

% % V1 Gray
% doOneFile(dirName, d.v1_c_gray(1), 10)
% set(gcf,'name','V1, control Gray 1')
% 
% doOneFile(dirName, d.v1_a_gray(1), 10)
% set(gcf,'name','V1, Atropine Gray 1')

%%
    function es_spk = doOneFile(dirName, x, chn, cell_list)
        if nargin<3
            chn = 20;
        end
        if nargin<4
            cell_list = 1:length(x.units);
        end
        fid = fopen([dirName filesep x.filename]);
        lfpData = fread(fid,'int16');
        fclose(fid);
        lfpData = reshape(lfpData, 32, []);
        [es, procSpec] = LFP_power_only(lfpData, chn, 1250, [2 1]);
        
        exptInfo.range_low = 1;
        exptInfo.range_high = 95;
        params.Fs = 1250;
        params.fpass=[exptInfo.range_low exptInfo.range_high];
        params.tapers=[5 9];
        clear spkTimes
        for icell = 1:length(x.units)
            spkTimes{icell} = x.units(icell).unit_spiketimes/30000;
            [C,phi,S12,SLFP,es_spk.powA(icell,:,:),es_spk.t,es_spk.freq]=...
                cohgramcpt(lfpData(20,:)',spkTimes{icell},[2 1],params);
        end
%         spkTimes = cell2struct(spkTimes);
%         [es_spk.powA, es_spk.t, es_spk.freq] = mtspecgrampt(spkTimes, [2 1], params);
        
        figure
        ax(1) = subplot(311);
        f = es.freq>25 & es.freq<75;
        imagesc(es.t, es.freq(f), log((ones(size(es.powA,1),1)*es.freq(f)).*es.powA(:,f))')
        xlims = xlim;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        ylabel('Freq (Hz)');
        title('LFP spectrogram');
        
        ax(2) = subplot(312);
        f = es_spk.freq>25 & es_spk.freq<75;
%         imagesc(es_spk.t, es_spk.freq(f), sq(nanmean((es_spk.powA(cell_list,:,f)),1))');
        imagesc(es_spk.t, es_spk.freq(f), sq(nanmean((es_spk.powA(:,:,f)),1))');
%         xlims = xlim;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        ylabel('Freq (Hz)');
        title('Spiking spectrogram');
        
        ax(3) = subplot(313);
        es_spk.smthBallSpd = smthInTime(x.locomotion.speed,90, 250);
        plot(x.locomotion.timestamps/30000,...
            smthInTime(x.locomotion.speed,90, 250))
        axis tight
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none','xlim',xlims);
        ylabel('Run speed (cm/s)')
        xlabel('Time (s)');
        title('Running speed')
        linkaxes(ax,'x');
    end
end