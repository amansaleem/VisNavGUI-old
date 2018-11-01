function [map,x] = fast1Ddelaymap(X,Y,dx,scaling,nbXbinsmth,visualspeed,Nquantile)
if ~isempty(X)
    if nargin < 3
        dx = 2.5;%(max(X)/nbbins);
        scaling = 1;
    end
    if nargin < 4
        scaling = 1;
    end       
    nbbins = round(max(floor(X/dx)+1));
    x = 0:dx:dx*(nbbins-1);
    if nargin < 5
        nbXbinsmth = nbbins;
    end  
    
    Xpos = floor(X/dx)+1;
    map = zeros(Nquantile,max(Xpos));
    ntrajbins = max(Xpos);
    spdquantilelim = zeros(ntrajbins,2);       
    for tt = 1:Nquantile
        for xx = 1:ntrajbins
%             spdquantilelim(xx,1) =  mean(visualspeed(Xpos == xx)) + 5*(tt-1-Nquantile/2);%quantile(visualspeed(Xpos == xx),(tt-1)/Nquantile);
%             spdquantilelim(xx,2) = mean(visualspeed(Xpos == xx)) + 5*(tt-Nquantile/2);%quantile(visualspeed(Xpos == xx),tt/Nquantile);
            spdquantilelim(xx,1) =  quantile(visualspeed(Xpos == xx),(tt-1)/Nquantile);
            spdquantilelim(xx,2) = quantile(visualspeed(Xpos == xx),tt/Nquantile);
            %                     spdquantilelim(xx,1) = quantile(obj.data.es.smthTrajSpd(obj.Subset.all & itraj == xx),(phs-1)/obj.Bayes.nthetaphsbins);
            %                     spdquantilelim(xx,2) = quantile(obj.data.es.smthTrajSpd(obj.Subset.all & itraj == xx),phs/obj.Bayes.nthetaphsbins);
        end
        spdidx = visualspeed>=spdquantilelim(Xpos,1) & visualspeed<spdquantilelim(Xpos,2);
        matI = zeros(numel(Xpos),max(Xpos));
        for xx = 1:max(Xpos)
            matI(Xpos == xx,xx) = 1;
        end
        matO = Y.*scaling;
        matH=matI'*matI;
        lambda = 0;
        eyeRegMat = eye(size(matH));
        matH = matH + eyeRegMat*lambda*trace(matH);
        maptemp = matH\(matI(spdidx,:)'*matO(spdidx));
        maptemp = special_smooth_1d([maptemp;maptemp;maptemp],1/(3*nbXbinsmth),0,(3*nbbins));
        map(tt,:) =  maptemp(nbbins+1:2*nbbins);%matH\(matI'*matO);
        map(tt,:) = map(tt,:);
    end
else
    map = [];
    x = [];
end

end