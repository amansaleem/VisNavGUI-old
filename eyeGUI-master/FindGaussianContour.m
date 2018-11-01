% fits multivariate gaussian to COM of pupil
% to find COM, first take pixel of max darkness, zoom in there
% compute COM of zoomed in region
% recenter box on new COM
% fit gaussian
function [params] = FindGaussianContour(r,tpt)

params.xy     = [];
params.area   = 0;
params.mu     = [NaN NaN];
params.isgood = 0;

frame = r.fr(:,:,tpt);
r.nX  = size(frame,1);
r.nY  = size(frame,2);

%
% zero out pixels < saturation level
fr    = frame;
fr    = 255-fr;
fr    = max(0, fr-(r.sats));

% fr(fr > 0) = 1;
CC = bwconncomp(double(fr>0),4);
maxsize = 0;
largestobj = 1;
delpixel = [];
for k = 1:numel(CC.PixelIdxList)
    [ix,iy] = ind2sub(size(fr),CC.PixelIdxList{k});
    if sum(ismember(ix,[1:4 size(fr,1)-3:size(fr,1)])) == 0 && sum(ismember(iy,[1:4 size(fr,1)-3:size(fr,1)])) == 0 
        if numel(CC.PixelIdxList{k}) > maxsize
            largestobj = k;
            maxsize = numel(CC.PixelIdxList{k});
        end
    end
end
for k = 1:numel(CC.PixelIdxList)
    if k ~= largestobj
        delpixel = [delpixel CC.PixelIdxList{k}(:)'];
    end
end
fr(delpixel) = 0;
[Y,X] = find(fr>0);
if numel(X) > 2
    try
        K = convhull(X,Y);
        p = pdist([X(K) Y(K)]);
        Z = squareform(p);
        [xpts,ypts] = find(Z == max(Z(:)));
        p1 = [X(K(xpts(1))) Y(K(xpts(1)))];p2 = [X(K(ypts(1))) Y(K(ypts(1)))];
        Xsym = zeros(size(X));
        Ysym = zeros(size(Y));
        m = (p1(2)-p2(2))/(p1(1)-p2(1));n = p1(2) - m*p1(1);
        for i = 1:numel(X)
            Xsym(i) = (X(i) + m*Y(i) - m*n)/(m^2 + 1);
            Ysym(i) = m*Xsym(i) + n;
            Xsym(i) = 2*Xsym(i) - X(i);
            Ysym(i) = 2*Ysym(i) - Y(i);
        end
        
        IND = 1:numel(fr);
        [I,J] = ind2sub(size(fr),IND);
        [in,on] = inpolygon(I,J,Y(K),X(K));
        [insym,onsym] = inpolygon(I,J,Ysym(K),Xsym(K));
        in = in | insym;
        on = on | onsym;
        fr = 0*fr;
        for i = 1:numel(I)
            if on(i) || in(i)
                fr(I(i),J(i)) = 1;
            end
        end
        w = gausswin(20);
        fr = conv2(fr,w*w','same');
        
        % fr(X(K),Y(K)) = 1;
        %fr(fr>40) = 255;
        %fr    = 255 - fr;
        %imagesc(fr)
        % find pixel of max brightness
        [~,ix] = max(fr(:));
        [ix,iy] = ind2sub(size(fr),ix);
        
        % find com in window of 1/3 ROI size
        ixinds = 1:r.nX;%ix + [-1*round(r.nX/2):round(r.nX/2)];
        ixinds(ixinds>r.nX | ixinds<1) = [];
        iyinds = 1:r.nY;% iy + [-1*round(r.nY/2):round(r.nY/2)];
        iyinds(iyinds>r.nY | iyinds<1) = [];
        iyinds = repmat(iyinds(:), 1, numel(ixinds));
        ixinds = repmat(ixinds(:)', size(iyinds,1), 1);
        %
        ix     = sub2ind(size(fr), ixinds, iyinds);
        com    = [sum(ixinds(:).*fr(ix(:))) sum(iyinds(:).*fr(ix(:)))] /sum(fr(ix(:)));
    catch
        com = NaN;
    end
    % recenter box on com
    if ~isnan(com(1))
        ix     = round(com(1));
        iy     = round(com(2));
        ixinds = ix + [-1*round(r.nX/4):round(r.nX/4)];
        ixinds(ixinds>r.nX | ixinds<1) = [];
        iyinds = iy + [-1*round(r.nY/4):round(r.nY/4)];
        iyinds(iyinds>r.nY | iyinds<1) = [];
        iyinds = repmat(iyinds(:), 1, numel(ixinds));
        ixinds = repmat(ixinds(:)', size(iyinds,1), 1);
        ix     = sub2ind(size(fr), ixinds, iyinds);
        
        if sum(fr(ix(:))>0) > 1
            params = FitMVGaus(iyinds(:), ixinds(:), fr(ix(:)), r.thres);
            params.isgood = 1;
            params.fr = fr;
        end
    end    
end
