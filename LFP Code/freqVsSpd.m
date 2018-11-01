function [output] = freqVsSpd(pow_es, ranges)

if nargin < 2
    ranges = [3 12]';
end

subset = pow_es.freq>ranges(1) & pow_es.freq<ranges(2);

output.spd = pow_es.spdBins;
output.conf= pow_es.timeInBins;

for ispd = 1:size(pow_es.powA,1)
    output.freqA(ispd) = sum(pow_es.powA(ispd,subset) .* pow_es.freq(subset))...
        ./sum(pow_es.powA(ispd,subset));
    output.freqB(ispd) = sum(pow_es.powB(ispd,subset) .* pow_es.freq(subset))...
        ./sum(pow_es.powB(ispd,subset));
    output.freqAB(ispd) = sum(pow_es.powAB(ispd,subset) .* pow_es.freq(subset))...
        ./sum(pow_es.powAB(ispd,subset));
    output.coh(ispd) = sum(pow_es.coherence(ispd,subset) .* pow_es.freq(subset))...
        ./sum(pow_es.coherence(ispd,subset));
    
    output.areaA(ispd) = sum(pow_es.powA(ispd,subset));
    output.areaB(ispd) = sum(pow_es.powB(ispd,subset));
    output.areaAB(ispd) = sum(pow_es.powAB(ispd,subset));
    output.areaCoh(ispd) = sum(pow_es.coherence(ispd,subset));
    if isfield(pow_es, 'booterror')
        for iter = 1:size(pow_es.booterror.powA,1)
            output.error.freqA(ispd,iter) = sum(reshape(pow_es.booterror.powA(iter, ispd,subset),1,[])...
                .* pow_es.freq(subset)) ./sum(reshape(pow_es.booterror.powA(iter, ispd,subset),1,[]));
            output.error.freqB(ispd,iter) = sum(reshape(pow_es.booterror.powB(iter, ispd,subset),1,[])...
                .* pow_es.freq(subset)) ./sum(reshape(pow_es.booterror.powB(iter, ispd,subset),1,[]));
            output.error.freqAB(ispd,iter) = sum(reshape(pow_es.booterror.powAB(iter, ispd,subset),1,[])...
                .* pow_es.freq(subset)) ./sum(reshape(pow_es.booterror.powAB(iter, ispd,subset),1,[]));
            output.error.coh(ispd,iter) = sum(reshape(pow_es.booterror.coherence(iter, ispd,subset),1,[])...
                .* pow_es.freq(subset)) ./sum(reshape(pow_es.booterror.coherence(iter, ispd,subset),1,[]));
        end
        output.std.freqA(ispd) =  std(output.error.freqA(ispd,:));
        output.std.freqB(ispd) =  std(output.error.freqB(ispd,:));
        output.std.freqAB(ispd) =  std(output.error.freqAB(ispd,:));
        output.std.coh(ispd) =  std(output.error.coh(ispd,:));
    end
end


% % To plot all the conditions
% plot(output_c0.spd(2:end), output_c0.freqA(2:end),'o-','Color',[.5 .5 .5],'linewidth',1)
% hold on;
% plot(output.spd(2:end), output.freqA(2:end),'o-','Color',[1 0 0],'linewidth',1)
% plot(output_OL.spd(2:end), output_OL.freqA(2:end),'o-','Color',[0 0 1],'linewidth',1)
% legend('Gray','Closed-loop','Open-loop')
% set(gca, 'color','none','TickDir','out','box','off', 'fontsize',14);