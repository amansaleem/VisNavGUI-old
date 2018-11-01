function [SInfo, SInfoperSpk, map, x] = SpatialInformation(X,Y,dx,scaling,nbXbinsmth,Fcircular)
if ~isempty(X)
    if nargin < 3
        dx = 2.5;%(max(X)/nbbins);
        scaling = 1;
    end
    if nargin < 4
        scaling = 1;
    end
    if nargin < 6
        Fcircular = true;
    end
    nbbins = max(floor(X/dx)+1);%round(max(floor(X/dx)+1));
    x = 1:dx:max(floor(X/dx)+1);
    if nargin < 5
        nbXbinsmth = [];
    end
    if isempty(nbXbinsmth)
        nbXbinsmth = nbbins;
    end
        
    Y = Y.*scaling;
    meanRate = mean(Y);
    scMap = full(sparse(floor(X/dx)+1, 1, Y', nbbins, 1));
    scMapsqr = full(sparse(floor(X/dx)+1, 1, (Y.^2)', nbbins, 1));
    occMap = full(sparse(floor(X/dx)+1, 1, 1, nbbins, 1));
    if Fcircular
        scMap = [scMap;scMap;scMap];
        scMapsqr = [scMapsqr;scMapsqr;scMapsqr];
        occMap = [occMap;occMap;occMap];
    end
    scMap = special_smooth_1d(scMap,1/(3*nbXbinsmth),0,(3*nbbins));
    occMap = special_smooth_1d(occMap,1/(3*nbXbinsmth),0,(3*nbbins));
    map = scMap./occMap;
    
    if Fcircular
        map = map(nbbins+1:2*nbbins);
        occMap = occMap(nbbins+1:2*nbbins);
    end
    
    SInfo = sum(occMap/sum(occMap).*map.*log2(map/meanRate))*max(round(X))/nbXbinsmth;
    SInfoperSpk = SInfo/meanRate;%SInfo/sum(occMap/sum(occMap).*map);
    
else
    map = [];
    SInfo = 0;
    SInfoperSpk = 0;
end

end