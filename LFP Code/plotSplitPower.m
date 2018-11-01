function plotSplitPower(P, P2, f, ranges)

if nargin<4
    ranges = [1 20; 20 47; 52 98]';
end

if isempty(P2);
    twoChn = 0;
else
    twoChn = 1;
end

numRanges = size(ranges,2);
numSpds   = size(P,1);

lws = 0.75:(3.25./numSpds):4;
lcs = 0:(0.75./numSpds):0.75;

% figure('Position',[500 450 500 500]);
for irange = 1:numRanges
    currRange = f>=ranges(1,irange) & f<=ranges(2,irange);
    
%     subplot(1,numRanges,irange)
    hold off;
    for iSpd = 1:numSpds
        plot(f(currRange), abs(mean(P(iSpd,currRange),1)),...
            'color',lcs(iSpd)*[1 1 1], 'linewidth',lws(iSpd));
        hold on;
        if twoChn
            plot(f(currRange), abs(mean(P2(iSpd,currRange),1)),...
                'color',lcs(iSpd)*[1 0 0], 'linewidth',lws(iSpd));
        end
    end
    set(gca, 'color','none','TickDir','out','box','off', 'fontsize',14);
%     set(gca,'YTick',[]);
%     set(gca,'YScale','log');
%     set(gca,'XScale','log');
    xlabel('Freq (Hz)');
    if irange==1
        ylabel('Power (a.u.)');
    end
    axis tight
    
%     subplot(2,numRanges,irange+numRanges)
%     hold off;
%     for iSpd = 1:numSpds
%         plot(f(currRange), abs(f(currRange).*mean(P(iSpd,currRange),1)),...
%             'color',lcs(iSpd)*[1 1 1], 'linewidth',lws(iSpd));
%         hold on;
%         if twoChn
%             plot(f(currRange), abs(f(currRange).*mean(P2(iSpd,currRange),1)),...
%                 'color',lcs(iSpd)*[1 0 0], 'linewidth',lws(iSpd));
%         end
%     end
%     set(gca, 'color','none','TickDir','out','box','off', 'fontsize',14);
%     xlabel('Freq (Hz)');
%     set(gca,'YTick',[]);
%     if irange==1
%         ylabel('Power (a.u.)');
%     end
%     axis tight
end