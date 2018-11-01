function predave = getCircularAverage(mat,amp_th0,maxtol,interpfactor)
    mat(isnan(mat)) = 0;
    
    if nargin < 4
        interpfactor = 1;
    end
    if interpfactor ~=1
        Prange = size(mat,1);
        Xrange = size(mat,2);
        mat = [mat(end,:);mat];
        mat = interp1(0:Prange,mat,interpfactor:interpfactor:Prange,'spline');
        if Xrange == 1
            mat = mat(:);
        end
    end
    
    Prange = size(mat,1);
    Xrange = size(mat,2);
    postfilt = zeros(Prange,1);
    postfilt(1:Prange) = (1:(Prange))'/(Prange)*(2*pi);
    
    x = 1:Prange;
    cutoff = 0;%1;
    [maxval, ~] = max(mat,[],1);
    [minval, ~] = min(mat,[],1);
    for i = 1:Xrange
        if amp_th0 ~= 0
            [~,igood,ampth_i] = findfield(mat(:,i),amp_th0);
            if ~isempty(igood)
                mat(~ismember(x,igood),i) = cutoff;%ampth_i;%
            end
        else
            ampth_i = 0;
        end
        mat(mat(:,i) < minval(i) + (1-maxtol)*(maxval(i)-minval(i)),i) = cutoff;%ampth_i;%
    end
    mat2 = mat;
%     mat2 = zeros(size(mat,1)+1,size(mat,2));
%     mat2(1:end-1,:) = mat;
%     mat2(end,:) = mat(1,:);
    a = sum((mat2'*cos(postfilt)),2);b = sum((mat2'*sin(postfilt)),2);
    predave = atan2(b,a);%predave = circ_mean(postfilt, mat);
    predave = mod(predave + 2*pi, 2*pi);
    predave = (interpfactor*(predave/(2*pi)*Prange));
    predave(sum(mat2,1)==0) = NaN;
end