function [map1, map2, map3, k1, k2, k3] = regressSpeed(es, subset, dangle)

if nargout>1
    proc_traj = 1;
else
    proc_traj = 0;
end
if nargin<3
    dangle = 45;
    
end
if subset
    t = es.smthBallSpd>1 & es.smthTrajSpd>1;
else
    t = true(size(es.smthBallSpd));
end

X1 =  es.smthBallSpd(t);
if proc_traj
    X2 = es.smthTrajSpd(t);
    X3 = cosd(dangle)*es.smthTrajSpd(t) + sind(dangle)*es.smthBallSpd(t);
end
t1 = ~isnan(X1) & ~isinf(X1);
if proc_traj
    t2 = ~isnan(X2) & ~isinf(X2);
    t3 = ~isnan(X3) & ~isinf(X3);
end

k1 = oneDimMap;
k1.kfold = 1;
k1.smth_win = 1;
k1.sampleRate = 1000;
if proc_traj
    k2 = oneDimMap;
    k2.kfold = 1;
    k2.smth_win = 1;
    k2.sampleRate = 1000;
    k3 = oneDimMap;
    k3.kfold = 1;
    k3.smth_win = 1;
    k3.sampleRate = 1000;
end

for ifreq = 1:length(es.freq)
    k1 = k1.trainSpikeMap(X1(t1), es.powB(t1,ifreq),1);
    map1(ifreq,:) = k1.model.tuning.respModel;
    
    if proc_traj
        k2 = k2.trainSpikeMap(X2(t2), es.powB(t2,ifreq),1);
        k3 = k3.trainSpikeMap(X3(t3), es.powB(t3,ifreq),1);
        map2(ifreq,:) = k2.model.tuning.respModel;
        map3(ifreq,:) = k3.model.tuning.respModel;
    end   
end