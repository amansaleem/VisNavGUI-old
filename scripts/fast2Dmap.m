function [map,x,y] = fast2Dmap(X,Y,Z,dx,dy,scaling,nbXsmthbins,nbYsmthbins,Fcircular,nbXbins,nbYbins)
if ~isempty(X)
    if nargin < 4
        dx = 2.5;%(max(X)/nbbins);
        scaling = 1;
    end
    if nargin < 5
        dy = 10;%(max(X)/nbbins);
        scaling = 1;
    end
    if nargin < 6
        scaling = 1;
    end 
    
    if nargin < 10
        nbXbins = round(max(floor(X/dx)+1));
        nbYbins = round(max(floor(Y/dy)+1));
    end
    x = 0:dx:dx*(nbXbins-1);
    y = 0:dy:dy*(nbYbins-1);
    
    if nargin < 7
        nbXsmthbins = nbXbins;
    end
    Z = Z.*scaling;
    
    scMap = full(sparse(floor(Y/dy)+1, floor(X/dx)+1, Z', nbYbins, nbXbins));
    occMap = full(sparse(floor(Y/dy)+1, floor(X/dx)+1, 1, nbYbins, nbXbins));
    if Fcircular
        scMap = repmat(scMap,[3 3]);%([scMap scMap scMap];
        occMap = repmat(occMap,[3 3]);%[occMap occMap occMap];
    end
    map = special_smooth_2d(scMap,[1/nbYsmthbins 1/nbXsmthbins],0,0,[nbYbins nbXbins])./special_smooth_2d(occMap,[1/nbYsmthbins 1/nbXsmthbins],0,0,[nbYbins nbXbins]);
    
    if Fcircular
        map = map(nbYbins+1:2*nbYbins,nbXbins+1:2*nbXbins);%map(:,nbXbins+1:2*nbXbins);
    end
else
    map = [];
    x = [];
    y = [];
end

