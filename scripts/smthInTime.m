function out = smthInTime(in, sampFreq, smthWin, type, subset, smthtype)
%
% Usage: out = smthInTime(in, sampFreq, smthWin, type, subset, smthtype)
% to smooth in time the signal 'in'
%
% sampFreq = the sampling freq
% smthWin  = the half width of the smoothing window (in ms), 0 gives no
% smoothing
% type     = Shape argument of conv (optional), 'same' is default
% subset   = logical of region to be considered same size as input, everything else comes out
% as NaNs
%
% Aman Saleem
% November 2013
%JUL: 06.11.2015 - temporary fix to enable median filtering: if smthWin is
%positive, we do gaussian filtering; if negative, we do median filtering.

if nargin<5 || isempty(subset)
    subset = true(size(in));
end
if nargin<4 || isempty(type)
    type = 'same';
end
if nargin<6 || isempty(smthtype)
    smthtype = 'gaussian';
end

if smthWin==0
    out = in;
elseif strcmp(smthtype,'gaussian')
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(smthWin/samp_int);
    if win == 0
        out = in;
    else
        sGrid = max([100 win*5]);
        s = (-sGrid:sGrid);

        sfilt = (1./(win*(sqrt(2*pi)))).*exp(-(s.^2)./(2*(win^2)));

        % padding the inputs
        pad = zeros(size(sfilt));
        pad_length = length(pad);
        if size(in,1)~=1
            flip_again = true;
            in = in';
            subset = subset';
        else
            flip_again = false;
        end
        pad_in = [pad in pad];
        pad_subset = [pad double(subset) pad];

        pad_out         = conv(pad_in, sfilt, type);
        norm_subset     = conv(pad_subset, sfilt, type);

        pad_out = pad_out./norm_subset;
        out = pad_out((pad_length+1) : (end-pad_length));
        out(~subset) = NaN;

        if flip_again
            out = out';
        end
    end
elseif strcmp(smthtype,'gaussian_1')
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(smthWin/samp_int);
    if win == 0
        out = in;
    else        
        sGrid = max([100 win*5]);
        s = (-sGrid:sGrid);

        sfilt = (1./(win*(sqrt(2*pi)))).*exp(-(s.^2)./(2*(win^2)));
        sfilt = sfilt/max(sfilt);
        % padding the inputs
        pad = zeros(size(sfilt));
        pad_length = length(pad);
        if size(in,1)~=1
            flip_again = true;
            in = in';
            subset = subset';
        else
            flip_again = false;
        end
        pad_in = [pad in pad];

        pad_out         = conv(pad_in, sfilt, type);

        out = pad_out((pad_length+1) : (end-pad_length));
        out(~subset) = NaN;

        if flip_again
            out = out';
        end
    end
elseif strcmp(smthtype,'median')
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(abs(smthWin)/samp_int);
    if win == 0
        out = in;
    else
        out = medfilt1(in,win);
    end
 elseif strcmp(smthtype,'boxcar')
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(abs(smthWin)/samp_int);
    if win == 0
        out = in;
    else
        out = conv(in,ones(1,win)/win,'same');%smooth(in,win);
    end
    if smthWin > 0
        out = circshift(out,floor(win/2));
        out(1:floor(win/2)) = out(floor(win/2)+1);
    elseif smthWin < 0
        out = circshift(out,-floor(win/2));
        out(end-floor(win/2):end) = out(end-(floor(win/2)+1));
    end
 elseif strcmp(smthtype,'boxcar_centered')
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(smthWin/samp_int);
    if win == 0
        out = in;
    else
        out = conv(in,ones(1,win)/win,'same');%smooth(in,win);  
    end
 elseif strcmp(smthtype,'boxcarsum')
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(abs(smthWin)/samp_int);
    if win == 0
        out = in;
    else
        out = conv(in,ones(1,win),'same');%smooth(in,win).* smooth(win*ones(size(in)),win);
    end
    if smthWin > 0
        out = circshift(out,floor(win/2));
        out(1:floor(win/2)) = out(floor(win/2)+1);
    elseif smthWin < 0
        out = circshift(out,-floor(win/2));
        out(end-floor(win/2):end) = out(end-(floor(win/2)+1));
    end
    elseif strcmp(smthtype,'boxcarsum_centered')
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(abs(smthWin)/samp_int);
    if win == 0
        out = in;
    else
        out = conv(in,ones(1,win),'same');%smooth(in,win).* smooth(win*ones(size(in)),win);
    end
end