% Input:  a structure with any subset of the following fields (a user will be prompted
%                                           for whatever details that are not supplied)
%         .animal - a string
%         .series - a string or a number
%         .selectedExperiments - a cell array of strings or a column vector of numbers
%         .selectedChannels - a row vector specifying which channels (according to the original .ns5 channel naming convention)
%                            to put into the .dat file;
%                            if negative numbers are provided, the variable is treated as listing channels to omit
%         .identifier - a string which is used in forming the output filename
%         .internalRefContact - the channel (according to the original .ns5 channel naming convention)
%                              that gets subtracted from all the other ones (optional, [] by default)
%
% ztitchns5s relies on scripts & functions in \\zserver\code\Spikes  and  \\zserver\code\CerebusTools
% 2012-04 MO
function ztitchns5s(inputStruct)


SetDefaultDirs; % this script is in \\zserver\code\Spikes\
global DIRS;

if nargin < 1 || isempty(inputStruct)
  inputStruct.dummyField = [];
end;
if ~isstruct(inputStruct)
  error('ztitchns5s receives a single argument which must be a structure');
end;


%-------
%Determine animal name
if ~isfield(inputStruct, 'animal') || isempty(inputStruct.animal)
  tmp = inputdlg('Please provide the animal name');
  if isempty(tmp)
    error('No animal name provided');
  end;
  inputStruct.animal = tmp{1};
end;

%-------
%Get a valid identifier
if ~isfield(inputStruct, 'identifier') || isempty(inputStruct.identifier)
  tmp = inputdlg('Please provide the identifier to use in the output filenames');
  if isempty(tmp)
    error('No identifier provided');
  end;
  inputStruct.identifier = tmp{1};
elseif ~ischar(inputStruct.identifier)
  error('No valid identifier provided');
end;

%-------
%Determine the series
if ~isfield(inputStruct, 'series') || isempty(inputStruct.series)
  d=dir([DIRS.Cerebus filesep inputStruct.animal]);
  if isempty(d)
    error(['Directory ' DIRS.Cerebus filesep inputStruct.animal ' not found']);
  end;
  d=d(cell2mat({d(:).isdir})); % throw out the files, leave the subdirectories
  d=d(3:end); % throw out . and ..
  if isempty(d)
    error(['Directory ' DIRS.Cerebus filesep inputStruct.animal ' has no subdirectories']);
  end;
  d = {d.name};
  [series,tmp] = listdlg('PromptString','Please select the series',...
    'SelectionMode','single',...
    'ListString',d);
  if tmp < 1
    error('Error with selecting series...');
  end;
  inputStruct.series = d{series}; % series was an index into the list of choices, now it's a string
  clear d series;
elseif isnumeric(inputStruct.series)
  inputStruct.series = num2str(inputStruct.series);
end;

%-------
%Determine the experiments
if ~isfield(inputStruct, 'selectedExperiments') || isempty(inputStruct.selectedExperiments)
  d=dir([DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series]);
  if isempty(d)
    error(['Directory ' DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series ' not found']);
  end;
  d=d(cell2mat({d(:).isdir})); % throw out the files, leave the subdirectories
  d=d(3:end); % throw out . and ..
  if isempty(d)
    error(['Directory ' DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series ' has no subdirectories']);
  end;
  d = {d.name};
  [selectedExperiments,tmp] = listdlg('PromptString',{'Please select the', 'experiments'},...
    'InitialValue',1:length(d),...
    'ListString',d);
  if tmp < 1
    error(['Error with selecting experiments...']);
  end;
  inputStruct.selectedExperiments = {d{selectedExperiments}}; % selectedExperiments were indexes into the list of choices
  clear d selectedExperiments;
elseif isvector(inputStruct.selectedExperiments)
  tmp = {};
  for i = 1:length(inputStruct.selectedExperiments)
    tmp{i} = num2str(inputStruct.selectedExperiments(i));
  end;
  inputStruct.selectedExperiments = tmp;
  clear tmp i;
end;

if ~isfield(inputStruct, 'medianfilter')
  tmp = inputdlg('Do you wanna apply a median filter across channels');
  if isempty(tmp)
    error('No valid parameter');
  end;
  inputStruct.medianfilter = str2num(tmp{1});  
end

if ~isfield(inputStruct, 'whitening')
  tmp = inputdlg('Do you wanna perform a spatial whitening of the data?');
  if isempty(tmp)
    error('No valid parameter');
  end;
  inputStruct.whitening = str2num(tmp{1});  
end

if inputStruct.medianfilter || inputStruct.whitening
    sampRate = 30000;
    lowFq = 500;
    highFq = 0.95 * .5 * sampRate;
  end

if ~isfield(inputStruct, 'splittetrodes')
  tmp = inputdlg('Do you wanna split in tetrodes');
  if isempty(tmp)
    error('No valid parameter');
  end;
  inputStruct.splittetrodes = str2num(tmp{1});  
end


%-------
% Verify that each experiment directory has at least one .ns5 file
% and build a list of all the .ns5 files that will be unified into the .dat file
ns5files2unify = [];
for e = 1:length(inputStruct.selectedExperiments)
  ns5files = dir([DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series filesep inputStruct.selectedExperiments{e} filesep '*.ns5']);
  if isempty(ns5files)
    error([DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series filesep inputStruct.selectedExperiments{e} filesep ' contains no ns5 file']);
  end;
  
  % There are 2 possibilities now: if the data was saved 'by experiment',
  % there is just one .ns5 file in the directory. If, however, the data
  % was saved 'by repeat' there will be many files, and we will order them
  % by their date&time.
  for j=1:length(ns5files)
    ns5files2unify(end+1).fullFilename = [DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series filesep ...
      inputStruct.selectedExperiments{e} filesep ns5files(j).name];
    ns5files2unify(end).expName = inputStruct.selectedExperiments{e};
    ns5files2unify(end).date = datenum(ns5files(j).date);
  end;
  clear ns5files datearray dateindex j;
end;
clear e;

%------- Check that all the ns5 files have exactly the same no. of channels
for i = 1:length(ns5files2unify)
    fileH = fopen(ns5files2unify(i).fullFilename);
    if fileH < 0
        error(['ztitchns5s: Failed to open ', ns5files2unify(i).fullFilename]);
    end;
    fseek(fileH, 28, 'bof');
    chN(i)=fread(fileH, 1, 'uint32');
    fclose(fileH);
end;
if max(chN) ~= min(chN)
    error('ztitchns5s: not all ns5 files have the same no. of channels, a condition that (for now) is not handled!');
end;
clear i fileH;

%------- Now we re-order all the ns5 files according to their date & time
[~,dateindex]=sort(cell2mat({ns5files2unify(:).date}));
% dateindex = sort(dateindex);
ns5files2unify = ns5files2unify(dateindex);
clear dateindex;
% We also reorganize inputStruct.selectedExperiments to reflect the chronological order
inputStruct.selectedExperiments = {};
inputStruct.selectedExperiments{1} = ns5files2unify(1).expName;
for i = 2:length(ns5files2unify)
  if ~strcmp(inputStruct.selectedExperiments{end}, ns5files2unify(i).expName)
    inputStruct.selectedExperiments{end+1} = ns5files2unify(i).expName;
  end;
end;

%------- Construct the name of the output .dat filename (which we are about to create) & check it does not exist already:
outputDir = [DIRS.multichanspikes filesep inputStruct.animal filesep inputStruct.series filesep];
if ~exist(outputDir, 'dir')
  mkdir(outputDir);
  disp(['Created directory ' outputDir]);
end;
% outputDatFilename = [outputDir filesep inputStruct.animal '_' inputStruct.series '_' inputStruct.identifier '.dat'];
outputDatFilename = [outputDir filesep inputStruct.animal '_s' inputStruct.series '_' inputStruct.identifier];
outputMatFilename = strcat(outputDatFilename, '.mat');
% if exist(outputDatFilename, 'file') || exist(outputMatFilename, 'file')
%   error(['.dat/.mat file with the name ' outputDatFilename(1:end-4) ' already exists?!']);
% end;

%------- Present the user the list of all the ns5 files that will be stitched together
disp(['The following files (in the order shown) are about to be stitched into ', outputDatFilename]);
disp('----------------------------------------------------------------------------------------------');
disp(char(ns5files2unify(:).fullFilename));
disp('==============================================================================================');

%------- Now, let's see what happens with the channels...
[arrayLayout, ~, ~] = MichiganGetLayout(inputStruct.animal,str2num(inputStruct.series),chN); 
% not sure what we are supposed to do with the other two outputs (Create a .probe file?)

if ~isfield(inputStruct, 'selectedChannels') || isempty(inputStruct.selectedChannels)
  [selectedChannels,tmp] = listdlg('PromptString',{'Please select the channels to', 'include in the output .dat file'},...
    'InitialValue',1:length(arrayLayout),...
    'ListString',num2str(arrayLayout(:)));
  if tmp < 1
    error('Error with selecting channels...');
  end;
  inputStruct.selectedChannels = arrayLayout(selectedChannels); % selectedChannels were indexes into the list of choices
  clear selectedChannels;
  % here selectedChannels are already ordered according to the geometry of the probe...
end;
if sum(inputStruct.selectedChannels > 0) > 0 && sum(inputStruct.selectedChannels < 0) > 0
  error('inputStruct.selectedChannels contains both positive and negative entries');
elseif sum(inputStruct.selectedChannels < 0) > 0 % inputStruct.selectedChannels is a list of channels to omit
  [~, i]=setdiff(arrayLayout(:), abs(inputStruct.selectedChannels)); % i are the indexes of the remaining channels
  inputStruct.selectedChannels = arrayLayout(sort(i)); % the channels in their geometric order with the channels to omit thrown out
  clear i;
else % inputStruct.selectedChannels are the channels the user wishes to leave
  [~, i]=setdiff(arrayLayout(:), setdiff(arrayLayout(:), inputStruct.selectedChannels));
  tmp = arrayLayout(sort(i)); % the channels in their geometric order
  tmp = reshape(tmp, [], 1);
  inputStruct.selectedChannels = [tmp; setdiff(inputStruct.selectedChannels, tmp)']; 
  % the above allows the user to include channels other than those belonging to the probe in the .dat
  clear i tmp;
end;
clear arrayLayout;
CHANNELS_ORDER = inputStruct.selectedChannels;

if ~isfield(inputStruct, 'internalRefContact')%the number is in physical channel on a probe, before sorting
  inputStruct.internalRefContact = [];
end;

%------- Now we're sort of finally ready to get to the real work of writing the .dat file
global pepNEV;
lims = []; % this variable will contain the number of samples in each ns5 file that goes into the output .dat
for f = 1:length(ns5files2unify)
  % open .ns5 file and load data
  [tmp,~,nchan] = nsopen(ns5files2unify(f).fullFilename); % this function is in \\zserver\code\CerebusTools\
  if tmp < 0
    error(['Error nsopening ', ns5files2unify(f).fullFilename]);
  end;
  clear data;
  for ichan = 1:nchan
    data(ichan,:) = int16(pepNEV.ns.Data.data(ichan,:));
  end;  
  if ~isempty(inputStruct.internalRefContact)
    for ichan = 1:nchan
      data(ichan,:) = data(ichan,:) - data(inputStruct.internalRefContact,:);
    end;
  end;
  data = data(inputStruct.selectedChannels,:); % this performs the reordering according to probe geometry
  if inputStruct.medianfilter || inputStruct.whitening     
      for ichan = 1:size(data,1)
          data(ichan,:) = LFPfilter(double(data(ichan,:)), lowFq, highFq, sampRate);
      end;
  end
  
  if inputStruct.medianfilter
      channelmedian = median(double(data),1);
      for ichan = 1:size(data,1)
        data(ichan,:) = int16(double(data(ichan,:)) - channelmedian);
      end
  end
  
  if inputStruct.whitening 
      noise_th = 2;
      stdch = std(double(data),[],2);
      stdchmat = repmat(stdch,1,size(data,2));
      idx = sum(data > -noise_th*stdchmat & data < noise_th*stdchmat,1);
      noiseidx = find(idx==32);
      
      NoiseCovMat = double(data(:,noiseidx))*double(data(:,noiseidx))'/size(data(:,noiseidx),2);
      if inputStruct.splittetrodes
          NoiseCovMattemp = zeros(size(NoiseCovMat));
          for tet = 1:floor(size(data,1)/4)
              NoiseCovMattemp((tet-1)*4+1:tet*4,(tet-1)*4+1:tet*4) = NoiseCovMat((tet-1)*4+1:tet*4,(tet-1)*4+1:tet*4);
          end
          NoiseCovMat = NoiseCovMattemp;
      end
      WhitenMat = (double(NoiseCovMat)\eye(size(NoiseCovMat)))^0.5;
      data = int16(WhitenMat*double(data).*repmat(stdch,1,size(data,2)));%int16(WhitenMat*double(data)*100);      
  end  
  
  lims(end+1) = int64(size(data,2)); %can be a pretty big number, not to loose precision with float representation, we make it int64
  if ~inputStruct.splittetrodes
      if f == 1
          outputFileH = fopen([outputDatFilename '.dat'], 'w');
      else
          outputFileH = fopen([outputDatFilename '.dat'], 'a');
      end
      if outputFileH < 0
          error(['Failed to open ', outputDatFilename ' for writing']);
      end;
      tmp=fwrite(outputFileH, data, 'int16');
      if tmp ~= size(data,1) * size(data, 2)
          error(['Problem writing into ', outputDatFilename]);
      end;
      fclose(outputFileH);
      tmp = dir([outputDatFilename '.dat']);
      if tmp.bytes ~= sum(lims)*length(inputStruct.selectedChannels)*2 % sanity check
          error('Some problem writing the .dat file: its actual size is incorrect!');
      end;
  else
      for tet = 1:floor(size(data,1)/4)
          if f == 1
              outputFileH = fopen([outputDatFilename '_tet' num2str(tet-1) '.dat'], 'w');
          else
              outputFileH = fopen([outputDatFilename '_tet' num2str(tet-1) '.dat'], 'a');
          end
          if outputFileH < 0
              error(['Failed to open ', outputDatFilename  '_tet' num2str(tet-1) ' for writing']);
          end;
          tmp=fwrite(outputFileH, data((tet-1)*4+1:tet*4,:), 'int16');
          if tmp ~= 4 * size(data, 2)
              error(['Problem writing into ', outputDatFilename]);
          end;
          fclose(outputFileH);
      end      
  end
  nevclose; % presumably/hopefully also closes the .ns5 (?)
end;


%------- We save some parameters about the stitching into an acompanying .mat file:
SELECTED_CHANNELS = inputStruct.selectedChannels;
INTERNAL_REF_CONTACT = []; % The number (according to the channel order in .dat) of the channel used for internal reference
if ~isempty(inputStruct.internalRefContact)
  INTERNAL_REF_CONTACT = find(inputStruct.selectedChannels == inputStruct.internalRefContact);
end;
SELECTED_EXPERIMENTS = [];
for i = 1:length(inputStruct.selectedExperiments)
  SELECTED_EXPERIMENTS(end+1) = str2num(inputStruct.selectedExperiments{i});
end;
save(outputMatFilename, 'lims', ...
  'CHANNELS_ORDER', 'INTERNAL_REF_CONTACT', ...
  'SELECTED_CHANNELS', 'SELECTED_EXPERIMENTS', 'ns5files2unify', 'inputStruct');

