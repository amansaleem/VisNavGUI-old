function [protocol, successflag] = ProtocolLoad( animal, iseries, iexp, screenflag)
% ProtocolLoad load protocol data structure
%
% protocol = ProtocolLoad( ExptTag ) loads the protocol associated with
% ExptTag.animal, ExptTag.iseries, ExptTag.iexp.
%
% protocol = ProtocolLoad( animal, iseries, iexp )
%
% protocol = ProtocolLoad( animal, iseries, iexp, 'loadscreen')
%  loads the screenlogs of the experiment to compute the exact tfs
%
% protocol = ProtocolLoad( animal, iseries, iexp, myscreen)
%  uses myscreen to compute the exact tfs
%
% [ protocol, successflag ] = ProtocolLoad without arguments, uses the
% variables chosen in the global PICK
%
% Example: protocol = ProtocolLoad( 'CATZ008', 5, 2 )
%
% Important note: protocol.seqnum has an interpretation that is not intuitive.
% e.g. protocol.seqnum(2,3) = 7 means that stim_2 was the 7th one to be 
% presented in rep_3.
%
% In other words: The row indicates the stimulus number, and the content
% indicates the sequence number in which it was shown. So if on repeat 2
% stimulus 12 (of 13) was shown as the 4th in the sequence, then
% seqnums(12,2) = 17.
%
% part of Spikes toolbox

% 2000-09 MC created
% 2001-10 VM added protocol.estfreqs and .estfreqs2, eliminated oldprotocol in the input, added screenflag
% 2001-12 MC streamlined (was buggy if screendir does not exist)
% 2002-04 VB made sure loading screen log returns to original directory
% 2002-08 VB added PICK support
% 2003-02 VM added field sfnyquist and screenflag now can be the screen structure itself
% 2003-02 VB an empty structure is now returned if no protocol file found
% 2003-03 VM made it fit for traces
% 2003-07 VB added new v files.
% 2003-08 VB corrected bug. returned error when converting sf in experiments with only blank stimuli.
% 2003-10 MC corrected bug. returned error when reading an xfile that it does not know about.
% 2004-02 VB load screens if used with no arguments
% 2004-11 VB commented out disp statements to make code faster
% 2005-01 AB change in where to look for xfiles
% 2006-05 MC first argument can be ExptTag with appropriate fields.
% 2006-06 VB issues warning if could not find screen log
% 2006-08 MC loads the saved Protocol.mat file if it exists
% 2007-09 MC deals with PRIMING type
% 2008-09 LB inserted code to check for missing .mat files (bong on zazzera)
% 2008-10 MC more thorough check of which files are present (including Michigan)
% 2009-02 AZ (UNDONE) Append '\Cat' or '\Mouse' to DIRS.data, based on first letter of animal
% 2009-05 MC moved code that tests for file existence to ProtocolAbsentFiles
% 2010-02 MC added field "notes" to deal with buggy data...

global DIRS
global PICK

successflag = 1;

oldprotocol = ProtocolInitialize;

if nargin == 1
    ExptTag = animal; % the first argument was ExptInfo
    animal  = ExptTag.animal;
    iseries = ExptTag.iseries;
    iexp 	= ExptTag.iexp;
end

if nargin < 1
    animal = PICK.animal;
    iseries = PICK.iseries;
    iexp = PICK.iexp;
    screenflag = 'loadscreen';
end

if nargin < 4 && nargin > 0
    screenflag = 'donotload';
end

if ~isfield(DIRS,'xfiles')
    DIRS.xfiles = '';
end

ProtocolFileName = fullfile(DIRS.data,animal,int2str(iseries),int2str(iexp), 'Protocol.mat');

if exist( ProtocolFileName, 'file' )
    
    foo = load( ProtocolFileName ); % loads Protocol
    if isfield(foo,'protocol')
        % one of the new protocols, made by mpep...
        protocol = foo.protocol;
        [nstim,nrepeats] = size(protocol.seqnums);
        seqnums = zeros(nstim,nrepeats); % the corrected sequence
        for istim = 1:nstim
            for irepeat = 1:nrepeats
                seqnums(istim,irepeat) = find(protocol.seqnums(:,irepeat)==istim) + nstim*(irepeat-1);
            end
        end
        protocol.seqnums = seqnums;
%         for repIdx = 2:size(protocol.seqnums,2)
%             protocol.seqnums(:,repIdx) = protocol.seqnums(:,repIdx) + (repIdx-1)*size(protocol.seqnums,1);
%         end
        protocol.nstim = []; % so the field exists...
    elseif isfield(foo,'P')
        protocol = foo.P;
    else
        protocol = foo.Protocol;
    end
    
     if isempty(protocol.nstim)
        protocol.nstim = protocol.npfilestimuli;
     end

else
    
    fprintf('File %s does not exist!!!\n',ProtocolFileName);
    
    % BEGINNING OF CODE THAT SHOULD IN PRINCIPLE NEVER RUN AFTER AUGUST 2006
    % it runs if the protocol file has not been saved...
    fprintf('ProtocolLoad: apparently the protocol file has not been saved\n');
    
    %  read the pfile

    pfilename = sprintf('%s_%d_%d.p',animal, iseries, iexp);
    fullpfilename1 = fullfile(DIRS.data,animal,int2str(iseries),int2str(iexp),pfilename);
    fullpfilename2 = fullfile(DIRS.data,animal,int2str(iseries),pfilename);

    if exist(fullpfilename1,'file')==2
        fullpfilename = fullpfilename1;
    elseif exist(fullpfilename2,'file')==2
        fullpfilename = fullpfilename2;
    else
        fullpfilename = [];
    end

    if isempty(fullpfilename)
        protocol = struct([]);
        successflag = 0;
        errordlg(['Could not find pfile ' fullpfilename1 ' or '  fullpfilename2],...
            'Bad news from Spikes', 'modal');
        return;
    end

    pfile = fopen(fullpfilename,'r');

    if pfile < 0
        protocol = struct;
        successflag = 0;
        errordlg(['pfile ' fullpfilename ' exists but I could not open it'],'Spikes');
        return;
    end

    % disp([ 'Reading ' fullpfilename ]);

    protocol.xfile = fscanf(pfile,'%s\n',1);
    % if there is a path in the protocol.xfile name, remove it
    pathend = findstr(protocol.xfile,'/');
    if ~isempty(pathend)
        protocol.xfile = protocol.xfile(pathend+1:end);
    end

    readflag = 1;
    protocol.adapt.flag = 0;
    protocol.priming.flag = 0;
    while readflag
        flag = fscanf(pfile,'%s\n',1);
        switch flag
            case 'ADAPT'
                disp(' (this is an adaptation experiment) ');
                protocol.adapt.flag = 1;
            case 'PRIMING'
                disp(' (this is a priming experiment) ');
                protocol.priming.flag = 1;
            otherwise
                readflag = 0;
        end
    end

    protocol.npfilestimuli 	= str2num(flag);	% n of stimuli. can be protocol.nstim+2 if adapt exp...
    if isempty(protocol.npfilestimuli), protocol.npfilestimuli = 0; end
    
    protocol.npars 	= fscanf(pfile,'%d ',1);
    whatisthis 		= fscanf(pfile,'%d',1);
    screendistance	= fscanf(pfile,'%f ',1); % obsolete
    whatisthis 		= fscanf(pfile,'%d ',1);
    whatisthis 		= fscanf(pfile,'%d',1);

    pars = fscanf(pfile,'%d',protocol.npfilestimuli*protocol.npars);

    %--- done reading the parameter file
    fclose(pfile);

    if size(pars)~= protocol.npfilestimuli*protocol.npars
        protocol = struct;
        successflag = 0;
        disp('Pfile seems corrupted');
        return;
    end

    protocol.pars = zeros( protocol.npars, protocol.npfilestimuli); % important because it defines the size
    protocol.pars(:) = pars;

    % if needed, read the xfile
    
    if isempty(oldprotocol.xfile) || ~strcmp(oldprotocol.xfile,protocol.xfile)

        x = XFileLoad( protocol.xfile );
        
        protocol.xfile = x.name;
        protocol.npars = x.npars;
        protocol.parnames = x.parnames;
        protocol.pardefs = x.pardefs;
        
        if isempty(x.name)
            successflag = 0;
            disp('Could not read the x file');
            return
        end
    end

    protocol.nrepeats = 0; % added by MC Sept 5 2006
    
% END OF CODE THAT SHOULD IN PRINCIPLE NEVER RUN AFTER AUGUST 2006

end

if protocol.adapt.flag
    protocol.nstim = protocol.npfilestimuli-2;
else
    protocol.nstim = protocol.npfilestimuli;
end

protocol.animal  = animal;
protocol.iseries = iseries;
protocol.iexp    = iexp;

%% Fields that are conditional on the x file

% define these to be blank
% (in some cases, they will be overwritten)
tf2      = [];
isf      = [];
sfconv   = [];
idiam    = [];
diamconv = [];
nif     = [];

switch protocol.xfile
    case 'snddeltas.x'
        protocol.blankpars = [4 5];
        tf10pars = [];
    case 'visdriftsin.x'
        protocol.blankpars = [4]; % contrast
        protocol.ncycles 		= NaN*ones(1,protocol.nstim);
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        % FOR THE ANIMALS WITH THE SPATIAL FREQUENCY BUG
        buggyanimals = {...
            'CAT001', 'CAT002', 'CAT003','CATZ004','CATZ005','CATZ006',...
            'CATZ007','CATZ008','CATZ009','CATZ010','CATZ011','CATZ012'};
        if any(strmatch(upper(animal),buggyanimals,'exact'))
            sf 	= protocol.pars(3,:);
            ori 	= protocol.pars(5,:);
            diam 	= protocol.pars(8,:);
            actualsf = FixSpatFreqBug(sf,ori,diam);
            if any(actualsf~=sf)
                disp('----------------> Spatial frequencies were corrected to reflect the spatial frequency bug');
                protocol.pars(3,:) = actualsf;
            end
        end
        nif = 1;
        isf = [3];
        sfconv = [10];
        idiam = [8];
        diamconv = [10];

    case 'visdriftsin100.x'

        % copied from visdriftsin.x
        protocol.blankpars = [4]; % contrast
        protocol.ncycles 		= NaN*ones(1,protocol.nstim);
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        % FOR THE ANIMALS WITH THE SPATIAL FREQUENCY BUG
        buggyanimals = {...
            'CAT001', 'CAT002', 'CAT003','CATZ004','CATZ005','CATZ006',...
            'CATZ007','CATZ008','CATZ009','CATZ010','CATZ011','CATZ012'};
        if any(strmatch(upper(animal),buggyanimals,'exact'))
            sf 	= protocol.pars(3,:);
            ori 	= protocol.pars(5,:);
            diam 	= protocol.pars(8,:);
            actualsf = FixSpatFreqBug(sf,ori,diam);
            if any(actualsf~=sf)
                disp('----------------> Spatial frequencies were corrected to reflect the spatial frequency bug');
                protocol.pars(3,:) = actualsf;
            end
        end
        nif = 1;
        isf = [3];
        sfconv = [100];
        idiam = [8];
        diamconv = [10];

    case 'visdriftsin2.x'
        protocol.blankpars = [4 11]; % contrast
        protocol.ncycles 		= NaN*ones(1,protocol.nstim);
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        nif = 2;
        tf2 = protocol.pars(9,:)/10;
        isf = [3 10];
        sfconv = [10 10];
        idiam = [8 15];
        diamconv = [10 10];

    case 'vis2luts2grats.x'
        protocol.blankpars = [5 14]; % contrast
        protocol.ncycles 		= NaN*ones(1,protocol.nstim);
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        nif = 2;
        tf2 = protocol.pars(11,:)/10;
        isf = [3 12];
        sfconv = [10 10];
        idiam = [9 10 18 19];
        diamconv = [10 10 10 10];

    case 'vis2luts2grats100.x'
        protocol.blankpars = [5 14]; % contrast
        protocol.ncycles 		= NaN*ones(1,protocol.nstim);
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        nif = 2;
        tf2 = protocol.pars(11,:)/10;
        isf = [3 12];
        sfconv = [100 100];
        idiam = [9 10 18 19];
        diamconv = [10 10 10 10];

    case 'visplaid.x'
        protocol.blankpars = [5 6]; % contrast
        protocol.pfilefreqs = protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        protocol.ncycles 		= floor(protocol.pfiledurs.*protocol.pfilefreqs);
        nif = 1;
        isf = [3 4];
        sfconv = [10 10];
        idiam = [11];
        diamconv = [10];

    case 'vis2grat2win.x'
        protocol.blankpars = [5 6]; % contrast
        protocol.pfilefreqs = protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs = protocol.pars(1,1:protocol.nstim)/10;
        protocol.ncycles 		= floor(protocol.pfiledurs.*protocol.pfilefreqs);
        nif = 1;
        isf = [3 4];
        sfconv = [10 10];
        idiam = [11 12 13 14];
        diamconv = [10 10 10 10];

    case 'vis2grat2win2.x'
        protocol.blankpars = [7 8]; % contrast
        protocol.pfilefreqs = protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs = protocol.pars(1,1:protocol.nstim)/10;
        protocol.ncycles 		= floor(protocol.pfiledurs.*protocol.pfilefreqs);
        nif = 1;
        isf = [3 4];
        sfconv = [10 10];
        idiam = [13 14 15 16];
        diamconv = [10 10 10 10];

    case 'visringach.x'
        protocol.blankpars = [];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        isf = [2];
        sfconv = [10];
        idiam = [6];
        diamconv = [10];

    case 'visringnsfnori.x'
        protocol.blankpars = [];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        isf = [9];
        sfconv = [10];
        idiam = [5];
        diamconv = [10];

    case 'visringnsfnori100.x'
        protocol.blankpars = [];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        isf = [9];
        sfconv = [100];
        idiam = [5];
        diamconv = [10];

    case 'visringlog.x'
        protocol.blankpars = [];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        isf = [10];
        sfconv = [100];
        idiam = [5];
        diamconv = [10];

    case 'visringlogad.x'
        protocol.blankpars = 2;
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;

    case 'viscatcam.x'
        disp('Warning!!! Positions x and y are screwed up!');
        protocol.blankpars = [];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= [];
        isf = [];
        sfconv = [];
        idiam = [];
        diamconv = [];

    case 'vismovie.x'
        protocol.blankpars = [];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= [];
        isf = [];
        sfconv = [];
        idiam = [];
        diamconv = [];

    case 'vismradapt.x'
        protocol.blankpars = [];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= [];
        nif = 1;
        isf = [3];
        sfconv = [10];
        idiam = [6];
        diamconv = [10];

    case 'vis2grat2win2flash.x'
        protocol.blankpars = [10 11];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= [];
        nif = 2;
        isf = [6 7];
        sfconv = [100 100];
        idiam = [16 17 18 19];
        diamconv = [10 10 10 10];

    case 'visgratgrat.x'
        protocol.blankpars = [6];
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;;
        protocol.ncycles 		= floor(protocol.pfiledurs.*protocol.pfilefreqs);
        nif = 1;
        isf = [7 10];
        sfconv = [100 100];
        idiam = [5];
        diamconv = [10];

    case 'visgratgratlut.x'
        protocol.blankpars = [6];
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;;
        protocol.ncycles 		= floor(protocol.pfiledurs.*protocol.pfilefreqs);
        nif = 1;
        isf = [7 10];
        sfconv = [100 100];
        idiam = [5];
        diamconv = [10];

    case 'visgratgauss.x'
        protocol.blankpars = [6 11];
        protocol.pfilefreqs 	= protocol.pars(5,1:protocol.nstim)/10;;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;;
        protocol.ncycles 		= floor(protocol.pfiledurs.*protocol.pfilefreqs);
        nif = 1;
        isf = [7];
        sfconv = [100];
        idiam = [4 10];
        diamconv = [10 10];

    case 'vissweep.x'
        protocol.blankpars = [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= [];
        protocol.ncycles 		= [];
        isf = [7];
        sfconv = [100];
        idiam = [13];
        diamconv = [10];

    case 'vissweep2grat.x'
        protocol.blankpars = [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= [];
        protocol.ncycles 		= [];
        isf = [14 15];
        sfconv = [100 100];
        idiam = [4];
        diamconv = [10];

    case 'vislumcongauss.x'
        protocol.blankpars = [7 8 11 12 13];
        protocol.pfilefreqs 	= protocol.pars(14,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        protocol.ncycles 		= floor(protocol.pfiledurs.*protocol.pfilefreqs);
        nif = 1;
        isf = [15];
        sfconv = [100];
        idiam = [4 5 9 10];
        diamconv = [10 10 10 10];

    case 'vdriftsin100.x'
        % copied from visdriftsin.x
        protocol.blankpars = [4]; % contrast
        protocol.ncycles 		= NaN*ones(1,protocol.nstim);
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        nif = 1;
        isf = [3];
        sfconv = [100];
        idiam = [8];
        diamconv = [10];

    case 'vringlog.x'
        protocol.blankpars = 2;
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        isf = [10];
        sfconv = [100];
        idiam = [5];
        diamconv = [10];

    case {'vringlogad.x', 'oglRingach.x'}
        protocol.blankpars = 2;
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        isf = [10];
        sfconv = [100];
        idiam = [5];
        diamconv = [10];

    case 'vsweep2grat.x'
        protocol.blankpars = [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= [];
        protocol.ncycles 		= [];
        isf = [14 15];
        sfconv = [100 100];
        idiam = [4];
        diamconv = [10];

    case 'v2luts1grat1tex.x'
        protocol.blankpars = [7 14]; % contrast of grating, outer diameter of noise
        protocol.ncycles 		= NaN*ones(1,protocol.nstim);
        protocol.pfilefreqs 	= protocol.pars(4,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        nif = 2;
        isf = [5];
        sfconv = [100];
        idiam = [9 10 13 14];
        diamconv = [10 10 10 10];

    case 'v2luts1grat1texn.x'
        protocol.blankpars = [7 14];  % contrast of grating, outer diameter of noise
        protocol.ncycles 		= NaN*ones(1,protocol.nstim);
        protocol.pfilefreqs 	= protocol.pars(4,1:protocol.nstim)/10;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        nif = 2;
         isf = [5];
        sfconv = [100];
        idiam = [9 10 13 14];
        diamconv = [10 10 10 10];

    case 'vismeisteradapt.x'
        protocol.blankpars = [];
        protocol.ncycles 		= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= [];
        isf = [10 11];
        sfconv = [100 100];
        idiam = [3];
        diamconv = [10];

    case 'vlutgrat.x'
        protocol.blankpars = 4;
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.ncycles 	= floor(protocol.pfiledurs.*protocol.pfilefreqs);

    case 'v2lut2gratdel.x'
        protocol.blankpars = [5 15];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.ncycles 	= floor(protocol.pfiledurs.*protocol.pfilefreqs);

    case 'v2lut2gratdelbar.x'
        protocol.blankpars = [5 15];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10;
        protocol.ncycles 	= floor(protocol.pfiledurs.*protocol.pfilefreqs);

    case 'vmovie2gratinter.x'

        protocol.blankpars = [6 20];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10; % ignoring 'tf2'
        protocol.ncycles 	= floor(protocol.pfiledurs.*protocol.pfilefreqs);

    case 'vmovie2gratinterok.x'

        protocol.blankpars = [6 20];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        protocol.pfilefreqs 	= protocol.pars(2,1:protocol.nstim)/10; % ignoring 'tf2'
        protocol.ncycles 	= floor(protocol.pfiledurs.*protocol.pfilefreqs);
    
        % modified 080603 LB    
    case {'oglMovie2Grat.x', 'oglMovie2GratLin.x'}
        protocol.blankpars = [6 20];
        protocol.pfiledurs = protocol.pars(1,1:protocol.nstim)/10;
        protocol.pfilefreqs 	= [protocol.pars(2,1:protocol.nstim)/10; protocol.pars(16,1:protocol.nstim)/10]; % ignoring 'tf2'
        protocol.ncycles 	= [floor(protocol.pfiledurs.*protocol.pfilefreqs(1,:)); floor(protocol.pfiledurs.*protocol.pfilefreqs(2,:))];
    
        % added by MC 2010-02-05 because drifting gratings in
        % oglMovie2GratLin were buggy!
        
        if strcmp(protocol.xfile,'oglMovie2GratLin.x') && ~all(all(protocol.pars([12 26],:)))
        fprintf(1,'--------------> %s %d %d HAD BUGGY STIMULI!!!\n',...
            protocol.animal,protocol.iseries,protocol.iexp);
        protocol.notes = 'IGNORE!!! ';
        protocol.pars = zeros(protocol.npars,protocol.nstim)+eps; % purposefully weird numbers
        end
    
    case 'vmovie2sequentialGrating.x'
 
        protocol.blankpars = [7 21];
        protocol.pfiledurs 	= protocol.pars(1,:)/10;
        protocol.pfilefreqs = protocol.pars(17,:)/10; % looking only at 'tf2'
        protocol.ncycles 	= [];
 
    % added 2007-09-11 by MC
   case 'vXYCMSeqGratVer2.x'
 
        protocol.blankpars = 13;
        protocol.pfiledurs 	= protocol.pars(1,:)/10;
        protocol.pfilefreqs = [];
        protocol.ncycles 	= [];

    % Added 070702 LB; LB added oglKalatsky.x 080222    
    case {'kalatsky.x', 'oglKalatsky.x'}
        
        protocol.blankpars = 4; % contrast
        protocol.pfiledurs = protocol.pars(1,:)/10; % duration
        protocol.pfilefreqs = protocol.pars(2,:)/1000; % temporal freq
        protocol.ncycles = [];
        
    % Added 2007-09-05 by RF
    case 'vrandpos.x'
        
        protocol.blankpars = 2; % contrast
        protocol.pfiledurs = protocol.pars(1,:)/10; % duration
        protocol.pfilefreqs = [];
        protocol.ncycles = [];
        % Added 2008-01-28 by LB   

    case {'flickDrift.x', 'oglFlickDrift.x', 'oglflickDrift.x', 'oglFlickDrift02.x'}
        % Added 071022 by LB; LB added oglFlickDrift.x 080222; LB added oglFlickDrift02.x 080423
        protocol.blankpars = 4; % contrast
        protocol.pfiledurs = protocol.pars(1,:)/10 + protocol.pars(8,:)/10; % duration motion + duration blank
        % 2007-10 MC modified (frequency is 1/duration):
        protocol.pfilefreqs = 1./(protocol.pars(6,:)/10 + protocol.pars(7,:)/10); % 1/(duration_on + duration_off)
        protocol.ncycles = [];

    case {'oglRandpos.x', 'oglRandpos_sparse.x','oglRandpos_sparse2.x','oglRandpos_sparse2_alpha.x'} % Added 2010-08 by ND (for oglRandpos_sparse)
        protocol.blankpars = 2; % contrast
        protocol.pfiledurs = protocol.pars(1,:)/10; % duration
        protocol.pfilefreqs = [];
        protocol.ncycles = [];
        
    case {'oglRandpos_dense.x','oglRandpos_dense2.x','oglRandpos_dense2_alpha.x'} % Added 2010-08 by ND (for oglRandpos_dense)
        protocol.blankpars = 2; % contrast
        protocol.pfiledurs = protocol.pars(1,:)/10; % duration
        protocol.pfilefreqs = [];
        protocol.ncycles = [];

    case 'oglTwoGratings.x' % added by ND 2010-09
        protocol.blankpars = [6 20];
        protocol.pfiledurs 	= protocol.pars(1,1:protocol.nstim)/10;
        protocol.pfilefreqs = protocol.pars(2,1:protocol.nstim)/10; % ignoring 'tf2'
        protocol.ncycles 	= floor(protocol.pfiledurs.*protocol.pfilefreqs);
        
    case 'stimWaveOutput.x' % added by ND 2011-07
        protocol.blankpars  = [2 5 7 9 11]; % [amp v1 v2 v3 v4] must be 0 for all five in order to be a blank
        protocol.pfiledurs  = protocol.pars(1,1:protocol.nstim)/10;
        protocol.pfilefreqs = protocol.pars(3,1:protocol.nstim);
        protocol.ncycles    = floor(protocol.pfiledurs.*protocol.pfilefreqs);

    case {'flickeringChecks.x'}
        protocol.blankpars = []; % protocol.blankpars=3 % changed by ND to treat contrast as the activepar
        protocol.pfiledurs = protocol.pars(4,:)/10;
        protocol.pfilefreqs = protocol.pars(2,:)/10;
        protocol.ncycles = [];
        
    case {'oglFlickeringChecks.x'} % added by ND 2011-03
        protocol.blankpars = []; % protocol.blankpars=3 % changed by ND to treat contrast as the activepar
        protocol.pfiledurs = protocol.pars(1,:)/10;
        protocol.pfilefreqs = protocol.pars(2,:)/10;
        protocol.ncycles = floor(protocol.pfiledurs.*protocol.pfilefreqs);

    case {'flickeringChecks_stepresp.x'}
        protocol.blankpars = []; % protocol.blankpars=2 % changed by ND to treat contrast as the activepar
        protocol.pfiledurs = protocol.pars(3,:)/10;
        protocol.pfilefreqs = 1./(protocol.pars(3,:)/10); % frequency is 1/duration
        protocol.ncycles = [];
 
    case {'stimColorBar.x'}
        protocol.blankpars = [10 11 12]; 
        protocol.pfiledurs = protocol.pars(1,:)/10;
        protocol.pfilefreqs = 1./(protocol.pars(1,:)/10); % frequency is 1/duration
        protocol.ncycles = [];
        
    case {'stimColorBar2.x'}
        protocol.blankpars = [11 12 13]; 
        protocol.pfiledurs = protocol.pars(1,:)/10;
        protocol.pfilefreqs = 1./(protocol.pars(1,:)/10); % frequency is 1/duration
        protocol.ncycles = [];
 
    case {'stimKalatsky_xpos.x', 'stimKalatsky_ypos.x'} %added by AP Jun 2011 
        protocol.blankpars = [9 10 11]; 
        protocol.pfiledurs = protocol.pars(1,:)/10;
%         protocol.pfilefreqs = 1./(protocol.pars(1,:)/10); % frequency is 1/duration
        protocol.ncycles = protocol.pars(1,:).*protocol.pars(2,:)/1000;
        
    case {'stimKalatsky2_xpos.x', 'stimKalatsky2_ypos.x'} %added by AP Jun 2011 
        protocol.blankpars = [6 7 8]; 
        protocol.pfiledurs = protocol.pars(1,:)/10;
%         protocol.pfilefreqs = 1./(protocol.pars(1,:)/10); % frequency is 1/duration
        protocol.ncycles = protocol.pars(1,:).*protocol.pars(2,:)/1000;   

    case {'stimKalatsky.x'} %added by AP Jun 2012 
        protocol.blankpars = [11 12 13]; 
        protocol.pfiledurs = protocol.pars(1,:)/10;
%         protocol.pfilefreqs = 1./(protocol.pars(1,:)/10); % frequency is 1/duration
        protocol.ncycles = protocol.pars(1,:).*protocol.pars(4,:)/1000; 
 
    case {'stimKalatskyWavesFlicker.x'} %added by AP Jun 2013
        protocol.blankpars = [11 12 13];
        protocol.pfiledurs = protocol.pars(1,:)/10;
        %         protocol.pfilefreqs = 1./(protocol.pars(1,:)/10); % frequency is 1/duration
        protocol.ncycles = protocol.pars(1,:).*protocol.pars(4,:)/1000;
%     case {'vmovie3sequentialGrating.x'}
%        protocol.blankpars = [15 21]; 
%         protocol.pfiledurs = protocol.pars(1,:)/10;
%         protocol.ncycles = [];
%          protocol.pfilefreqs = [] % frequency is 1/duration


    case {'vmovie3sequentialGrating.x'}
        
        protocol.blankpars = 21; % this assumes that if c2 is zero it is a blank...
        
    case {'stimFlashedBarWaveOutput.x'} %added by ND Apr 2012
        protocol.blankpars = [12 14 17];
        protocol.pfiledurs = protocol.pars(1,:)/10;
        protocol.pfilefreqs = (protocol.pars(14,:)~=0).*(protocol.pars(15,:)); % frequency of optical stimulus
        
    case {'stimWavePulseVisPulse.x'} %added by ND Apr 2013
        protocol.blankpars = [6 14 16];
        protocol.pfiledurs = protocol.pars(1,:)/10;
%         protocol.pfilefreqs = (protocol.pars(14,:)~=0).*(protocol.pars(15,:)); % frequency of optical stimulus
        protocol.pfilefreqs = (protocol.pars(6,:)~=0).*(protocol.pars(9,:)); % frequency of visual stimulus

    case {'stimSweepGrating.x'} %added by AB March 2013
        protocol.blankpars  = 5; % contrast
        protocol.pfilefreqs = 1./(protocol.pars(1,:)/10); % frequency is 1/duration
        protocol.blankstims = protocol.nstim;
        
    otherwise
        disp(['Loading unknown xfile type ' protocol.xfile]);
        protocol.blankpars   = [];
        protocol.ncycles 	= [];
        protocol.pfilefreqs 	= [];
        protocol.pfiledurs 	= protocol.pars(1,:)/10;
        
end

%% find the effective temporal frequencies


%----------- default values:
protocol.estfreqs = [];
protocol.estfreqs2 = [];

%----------- now let's try to do better than that:

if isstr(screenflag) & ((strcmpi(animal,'catz') & str2num(animal(end-2:end)) > 12) | ~strcmp(lower(animal),'catz')) ...
        & strcmp(screenflag,'loadscreen')
    if ~isempty(nif) & ~isempty(protocol.pfilefreqs)
        [motherdir,datadir] = fileparts(DIRS.data);
        logname = sprintf('%s/%s_%i_%i.mat',upper(animal),lower(animal),iseries,iexp);
        screendir = fullfile(motherdir,'screen logs',upper(animal));
        if exist(screendir,'dir')
            previousdir=pwd;
            cd(screendir);
            screenfile = fullfile(motherdir,'screen logs',logname);
            if exist(screenfile,'file')
                fprintf(1,'Loaded screen log file %s in directory %s\n',logname,screendir);
                load(screenfile);
                FrameRate = myscreen.FrameRate;
                protocol.estfreqs = FrameRate ./ round(FrameRate ./ (protocol.pfilefreqs*nif)) / nif;
                if ~isempty(tf2)
                    protocol.estfreqs2 = FrameRate ./ round(FrameRate ./ (tf2*nif)) / nif;
                end
            end
            cd(previousdir);
        else
            fprintf('WARNING: Could not find screen log file %s and/or directory %s\n',logname,screendir);
        end
    end

elseif isstruct(screenflag)
    FrameRate = screenflag.FrameRate;
    protocol.estfreqs = FrameRate ./ round(FrameRate ./ (protocol.pfilefreqs*nif)) / nif;
    if ~isempty(tf2)
        protocol.estfreqs2 = FrameRate ./ round(FrameRate ./ (tf2*nif)) / nif;
    end
    disp('Computed estfreqs from the screen logs in the input');
end

%% find the blank stimuli

if isempty(protocol.blankpars)
    protocol.blankstims = [];
else
    protocol.blankstims = ones(1,protocol.nstim); % start by assuming they are all blank...
    for ipar = protocol.blankpars
        protocol.blankstims = protocol.blankstims & ~protocol.pars(ipar,1:protocol.nstim); % dismiss the non blanks...
    end
    protocol.blankstims = find(protocol.blankstims);
end

%% find the largest 'nyquist' frequency

notblank = setdiff(1:protocol.nstim,protocol.blankstims);
nyquist1 = [];
nyquist2 = [];
% Convert all sf into cycles/deg and find the max
if ~isempty(isf) & ~isempty(notblank)
    [junk,sfconvbig] = meshgrid(notblank,sfconv);
    sfsameunits = protocol.pars(isf,notblank)./sfconvbig;
    nyquist1 = max(sfsameunits(:));
end
% Convert all diams in deg and find the min
if ~isempty(idiam) && ~isempty(notblank)
    [junk,diamconvbig] = meshgrid(notblank,diamconv);
    diamsameunits = protocol.pars(idiam,notblank)./diamconvbig;
    notzero = find(diamsameunits ~= 0);
    nyquist2 = (min(diamsameunits(notzero)))^-1;
end
protocol.sfnyquist = max(nyquist1,nyquist2);

%% find the active variables

nonblanks = setdiff(1:protocol.nstim,protocol.blankstims);
protocol.activepars = {};
for ipar = 1:protocol.npars
    ipars = protocol.pars(ipar,nonblanks);
    ndiffipars = length(unique(ipars));
    if ndiffipars>1
        iactivepar = 1;
        foundflag = 0;
        while ~foundflag & iactivepar<=length(protocol.activepars)
            jpar = protocol.activepars{iactivepar}(1);
            jpars = protocol.pars(jpar,nonblanks);
            if (	length(unique(jpars)) == ndiffipars & ...
                    length(unique([ipars;jpars]','rows')) == ndiffipars )
                foundflag = 1;
                protocol.activepars{iactivepar} = [ protocol.activepars{iactivepar} ipar ];
            else
                iactivepar = iactivepar+1;
            end
        end
        if ~foundflag
            protocol.activepars{length(protocol.activepars)+1} = ipar;
        end
    end
end

%% figure out a description

protocol.description = [];

if isfield(protocol,'notes')
    protocol.description = protocol.notes;
end

for iactivepar = 1:length(protocol.activepars)
    parlist = protocol.activepars{iactivepar};
    protocol.description = [ protocol.description num2str(length(unique(protocol.pars(parlist(1),:)))) ];

    if ~ischar(protocol.parnames{parlist(1)})
        warning('Something buggy in protocol.parnames');
        break
    end
    
    protocol.description = [ protocol.description ' ' protocol.parnames{parlist(1)} ];
    for par = parlist(2:end)
        protocol.description = [ protocol.description ',' protocol.parnames{par} ];
    end
    protocol.description = [ protocol.description ' ' ];
end
if protocol.adapt.flag
    protocol.description = [protocol.description ' ADAPT '];
end
protocol.description = [protocol.description '(' num2str(protocol.npfilestimuli) ', ' protocol.xfile ')'];

%% add a quick tag (MC 2010-03)

protocol.tag = sprintf('%s-%d-%d',protocol.animal,protocol.iseries,protocol.iexp);

%% All done

return

%% function FixSpatFreqBug

function actualsf = FixSpatFreqBug(sf,ori,diam)
% FixSpatFreqBug
%
% actualsf = FixSpatFreqBug(sf,ori,diam)
% ori is between 0 and 360
%
% example:
%
% sf 	= [ 2 2 2 2 2 2 2 2 ]; % cycles/deg
% diam 	= [ 8 8 8 8 8 8 8 8 ];
% ori 	= [ 0 30 45 60 90 180 270 360];
%
% actualsf = FixSpatFreqBug(sf,ori,diam)
%
% 2001-03 MC

f = sf/10; % cycles/deg
d = diam/10; % deg
theta = ori * pi/180; % direction of motion, in radians

% actualx = d.* min(1,abs(cos(theta))./abs(sin(theta)));
% actualy = d.* min(1,abs(sin(theta))./abs(cos(theta)));

if abs(cos(theta))<abs(sin(theta))
    actualx = d.* abs(cos(theta))./abs(sin(theta));
else
    actualx = d;
end

if abs(cos(theta))>abs(sin(theta))
    actualy = d.* abs(sin(theta))./abs(cos(theta));
else
    actualy = d;
end

maxperiod = sqrt(actualx.^2+actualy.^2);

actualperiod = min(maxperiod,1./f);

actualf = 1./actualperiod;

actualsf = round(10*actualf);
