function [output_es, pow_es] = powAnalysis(animal, iseries, expt_list, cond_tag, chns, savedata)

if nargin< 5
    savedata = false;
end
if nargin< 4
    cond_tag = 'none';
end

if nargin< 5
    chns = [20 40];
end

SetDirs;

es = VRLoadMultipleExpts(animal,iseries,expt_list, 'BEHAV_LFP', chns);

allCont = unique(es.contrast(~isnan(es.contrast)));
if length(allCont)==4
    cont = allCont(3);
elseif length(allCont)>1
    cont = allCont(end-1);
else
    cont = allCont;
end

switch cond_tag
    case 'none'
        t = ones(size(es.traj));
    case 'active'
        t = es.traj~=0 & es.contrast ~=0 & es.outcome==2;% & es.smthBallSpd>2;
    case 'notactive'
        t = es.traj~=0 & es.contrast ~=0 & es.outcome~=2;% & es.smthBallSpd>2;
    case 'gray'
        t = es.traj==0 | es.contrast ==0;
    case 'OL'
        t = es.iexp>200 & es.traj~=0 & es.contrast ~=0;
    case 'passive'
        t = es.traj~=0 & es.contrast ~=0 & (es.outcome==0 | es.outcome==1);
    case 'CL'
        t = es.traj~=0 & es.contrast ~=0 & es.iexp<200 & es.iexp>100;
    case 'lowContrast'
        t = es.traj~=0 & es.contrast < cont;% & es.smthBallSpd>2;
    case 'highContrast'
        t = es.traj~=0 & es.contrast > cont;% & es.smthBallSpd>2;
    case 'lowGain'
        t = es.traj~=0 & es.contrast > 0 & es.gain < 1;% & es.smthBallSpd>2;
    case 'highGain'
        t = es.traj~=0  & es.contrast > 0 & es.gain > 1;% & es.smthBallSpd>2;
    case 'lowRL'
        t = es.traj~=0 & es.contrast > 0 & es.roomLength < 1;% & es.smthBallSpd>2;
    case 'highRL'
        t = es.traj~=0  & es.contrast > 0 & es.roomLength > 1;% & es.smthBallSpd>2;
    
end

pow_es = getPowSpd(es, t);
output_es = freqVsSpd(pow_es);

if savedata
    save([animal '_s' num2str(iseries) '_' cond_tag],'pow_es','output_es','es','t','-v7.3')
end