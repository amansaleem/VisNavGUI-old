function [q_idx,fieldix,amp_th_i] = findfield(tuning1D,qthreshold,imax,Fquantile)
if nargin < 3
    imax = [];
end
if nargin < 4
    Fquantile = true;
end

if min(tuning1D) ~= max(tuning1D) && sum(isnan(tuning1D)) == 0
    if Fquantile
        q_idx = [];
        if qthreshold > 0
            if qthreshold > 1
                factor = qthreshold;
                qthreshold = 1;
            else
                factor = 1;
            end
            if qthreshold < 1
                q_th = quantile(tuning1D,qthreshold);
                q_idx = find(tuning1D <= q_th);
            else
                q_idx = 1:numel(tuning1D);
            end
            
            if numel(q_idx) > 1
                Xrange = numel(tuning1D);
                amp_th_i = factor*nanmean(tuning1D(q_idx));
                if isempty(imax)
                    [~,imax] = max(tuning1D);
                end
                if ~isnan(amp_th_i)
                    iup = find(diff(sign([tuning1D(:);tuning1D(:);tuning1D(:)]-amp_th_i))>0)+1;
                    idown = find(diff(sign([tuning1D(:);tuning1D(:);tuning1D(:)]-amp_th_i))<0);
                    if ~isempty(iup) && ~isempty(idown)
                        iup = iup(find(iup<=imax+Xrange,1,'last'));
                        idown = idown(find(idown>=imax+Xrange,1,'first'));
                        fieldix = mod((iup:idown)-1, Xrange)+1;
                    end
                end
            else
                q_idx = [];
                amp_th_i = min(tuning1D);
                fieldix = 1:numel(tuning1D);
            end
        else
            q_idx = [];
            amp_th_i = min(tuning1D);
            fieldix = 1:numel(tuning1D);
        end
    else
        q_idx = [];
        Xrange = numel(tuning1D);
        amp_th_i = qthreshold;
        if ~isnan(amp_th_i)
            iup = find(diff(sign([tuning1D(:);tuning1D(:);tuning1D(:)]-amp_th_i))>0)+1;
            idown = find(diff(sign([tuning1D(:);tuning1D(:);tuning1D(:)]-amp_th_i))<0);
            if ~isempty(iup) && ~isempty(idown)
                iup = iup(find(iup<=imax+Xrange,1,'last'));
                idown = idown(find(idown>=imax+Xrange,1,'first'));
                fieldix = mod((iup:idown)-1, Xrange)+1;
            end
        end
        q_idx = find(~ismember(1:numel(tuning1D),fieldix));
    end
else
    q_idx = 1:numel(tuning1D);
    amp_th_i = 0;
    fieldix = [];
end
end