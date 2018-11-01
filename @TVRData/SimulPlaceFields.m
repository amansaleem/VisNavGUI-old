function obj = SimulPlaceFields(obj,basedelay,modeltype,VistoDistfactor)

samplingrate = 60;
baselinegain = min(obj.data.es.gain);%mode(obj.data.es.gain);
cbase = find(obj.SubsetVal.contrast == mode(obj.data.es.contrast));
gbase = find(obj.SubsetVal.gain == mode(obj.data.es.gain));
rbase = find(obj.SubsetVal.roomlength == mode(obj.data.es.roomLength));
obase = find(obj.SubsetVal.outcome == 2);

maxdist = max(floor(obj.data.es.distTrav)+1);
maxtraj = max(floor(obj.data.es.traj)+1);

obj.data.es.spikeTrain = 0 * obj.data.es.spikeTrain;
for icell = 1:size(obj.data.es.spikeTrain,2)
    delay = round(basedelay/(1000/samplingrate));% + randi(round(basedelay)) - round(basedelay/2); 
%     placefieldmed = obj.maps1d{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel;
%     placefieldlow = circshift(obj.maps1d{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel,[0 -10]);
%     placefieldhigh = circshift(obj.maps1d{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel,[0 10]);
    if strcmp(modeltype,'Vision')
        lambdaplace =  obj.maps1d{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel(floor(obj.data.es.traj)+1)/samplingrate;
    elseif strcmp(modeltype,'Global Distance')
        baselinegain = min(obj.data.es.gain);
        lambdaplace =  obj.maps1d{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel(min(floor(obj.data.es.distTrav*baselinegain)+1,maxtraj))/samplingrate;
    elseif strcmp(modeltype,'Local Distance')
        baselinegain = mode(obj.data.es.gain);
        lambdaplace = zeros(size(obj.data.es.traj'));
        for xbins = 1:VistoDistfactor            
            Xidx = find(floor(obj.data.es.traj)+1 >= (xbins-1)*maxtraj/VistoDistfactor & floor(obj.data.es.traj)+1 < xbins*maxtraj/VistoDistfactor);
            lambdaplace(Xidx) =  obj.maps1d{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel(min(floor((xbins-1)*maxtraj/VistoDistfactor + (obj.data.es.traj(Xidx)-(xbins-1)*maxtraj/VistoDistfactor)./(obj.data.es.gain(Xidx)/baselinegain))+1,maxtraj))/samplingrate;
        end
    elseif strcmp(modeltype,'Mixt')
        lambdaplace =  VistoDistfactor*obj.maps1d{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel(floor(obj.data.es.traj)+1)/samplingrate +...
                       (1-VistoDistfactor)*obj.maps1d{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel(floor(obj.data.es.distTrav*baselinegain)+1)/samplingrate;
    end
    lambdaplace = circshift(lambdaplace',delay);
    
%     lambdaplace(obj.data.es.gain==0.4) = placefieldlow(floor(obj.data.es.traj(obj.data.es.gain==0.4))+1)/samplingrate;
%     lambdaplace(obj.data.es.gain==0.6) = placefieldhigh(floor(obj.data.es.traj(obj.data.es.gain==0.6))+1)/samplingrate;
    
    lambdaval = unique(lambdaplace);
    for tt=1:numel(lambdaval)
        ttidx = find(lambdaplace == lambdaval(tt));
        obj.data.es.spikeTrain(ttidx,icell) = poissrnd(lambdaval(tt),[numel(ttidx),1]);%randi(2, sum(ismember(round(obj.data.es.traj), xfield)), 1) - 1;
    end
end
end