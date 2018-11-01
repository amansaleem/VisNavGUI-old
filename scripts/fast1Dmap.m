function [map,x,mapstd,mapsem] = fast1Dmap(X,Y,dx,scaling,nbXbinsmth,Fcircular)
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
    nbbins = round(max(floor(X/dx)+1));%round(100/dx);%round(max(floor(X/dx)+1));
    x = 0:dx:dx*(nbbins-1);
    if nargin < 5
        nbXbinsmth = [];
    end
    if isempty(nbXbinsmth)
        nbXbinsmth = nbbins;
    end
        
    Y = Y.*scaling;
    scMap = full(sparse(floor(X/dx)+1, 1, Y', nbbins, 1));
    scMapsqr = full(sparse(floor(X/dx)+1, 1, (Y.^2)', nbbins, 1));
    occMap = full(sparse(floor(X/dx)+1, 1, 1, nbbins, 1));
    if Fcircular
        scMap = [scMap;scMap;scMap];
        scMapsqr = [scMapsqr;scMapsqr;scMapsqr];
        occMap = [occMap;occMap;occMap];
        map = special_smooth_1d(scMap,1/(3*nbXbinsmth),0,(3*nbbins))./special_smooth_1d(occMap,1/(3*nbXbinsmth),0,(3*nbbins));
        mapsqr =  special_smooth_1d(scMapsqr,1/(3*nbXbinsmth),0,(3*nbbins))./special_smooth_1d(occMap,1/(3*nbXbinsmth),0,(3*nbbins));
    else
        map = special_smooth_1d(scMap,1/(nbXbinsmth),0,(nbbins))./special_smooth_1d(occMap,1/(nbXbinsmth),0,nbbins);
        mapsqr =  special_smooth_1d(scMapsqr,1/(nbXbinsmth),0,(nbbins))./special_smooth_1d(occMap,1/(nbXbinsmth),0,(nbbins));
    end
%     map = scMap ./ occMap;
%     map(isnan(map)) = 0;
%     map = special_smooth_1d(map,1/nbbins,0,nbbins);
    
%     map = special_smooth_1d(scMap,1/nbbins,0,nbbins)./special_smooth_1d(occMap,1/nbbins,0,nbbins);
%     map = scMap./occMap;
    
    mapstd = (mapsqr - map.^2).^0.5;%std
    mapsem = ((mapsqr - map.^2).^0.5)./special_smooth_1d(occMap,1/(3*nbXbinsmth),0,(3*nbbins)).^0.5;%sem
   
%     lambdasmooth = 4;
%     map = smooth1D(scMap,lambdasmooth)./smooth1D(occMap,lambdasmooth);

% Z = [X Y/60];
% [F, c1, c2, con, H] = smoothhist2D_corrected(Z, lambdasmooth, [nbbins max(Y/60)], 1:nbbins, 1:max(Y/60), true, false);
% map = F;

%     map = map * scaling;
    if Fcircular
        map = map(nbbins+1:2*nbbins);
        mapstd = mapstd(nbbins+1:2*nbbins);
        mapsem = mapsem(nbbins+1:2*nbbins);
    end
%     Y0 = Y;
%     mapfull = zeros(max(Y0)+1,nbbins);
%     for i = 0:max(Y0)
%         Y = double(round(Y0) == i);
%     
%         scMap = full(sparse(floor(X/dx)+1, 1, Y', nbbins, 1));
%         occMap = full(sparse(floor(X/dx)+1, 1, 1, nbbins, 1));
%     
%         scMap = [scMap;scMap;scMap];
%         occMap = [occMap;occMap;occMap];
%         map = special_smooth_1d(scMap,1/(3*nbXbinsmth),0,(3*nbbins))./special_smooth_1d(occMap,1/(3*nbXbinsmth),0,(3*nbbins)); 
%         map = map(nbbins+1:2*nbbins);
%         
%         mapfull(i+1,:) = map;
%     end
else
    map = [];
    x = [];
    mapstd = [];
    mapsem = [];
end

end
