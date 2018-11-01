function run_expt_list

ianim = 1;

% animal(ianim).name = 'M130405a_BALL';
% animal(ianim).series = 801;
% animal(ianim).chns = [10 40];
% 
% animal(ianim).tags(1).tag = 'CL';
% animal(ianim).tags(1).expt_list = 106:107;
% 
% animal(ianim).tags(2).tag = 'OL';
% animal(ianim).tags(2).expt_list = 201:203;
% 
% animal(ianim).tags(3).tag = 'gray';
% animal(ianim).tags(3).expt_list = 106:107;
% 
% ianim = ianim +1;
% %
% animal(ianim).name = 'M130702_BALL';
% animal(ianim).series = 807;
% animal(ianim).chns = [10 40];
% 
% animal(ianim).tags(1).tag = 'active';
% animal(ianim).tags(1).expt_list = 103:107;
% 
% animal(ianim).tags(2).tag = 'OL';
% animal(ianim).tags(2).expt_list = 201:203;
% 
% animal(ianim).tags(3).tag = 'gray';
% animal(ianim).tags(3).expt_list = 112;
% 
% animal(ianim).tags(4).tag = 'passive';
% animal(ianim).tags(4).expt_list = 103:107;
% 
% ianim = ianim +1;
%%
animal(ianim).name = 'M130703_BALL';
animal(ianim).series = 809;
animal(ianim).chns = [10 40];

animal(ianim).tags(1).tag = 'CL';
animal(ianim).tags(1).expt_list = 106:109;

animal(ianim).tags(2).tag = 'OL';
animal(ianim).tags(2).expt_list = 201:202;

animal(ianim).tags(3).tag = 'gray';
animal(ianim).tags(3).expt_list = [106 113];

animal(ianim).tags(4).tag = 'passive';
animal(ianim).tags(4).expt_list = 103:107;

animal(ianim).tags(5).tag = 'active';
animal(ianim).tags(5).expt_list = 103:107;

ianim = ianim +1;
%%
animal(ianim).name = 'M130918_BALL';
animal(ianim).series = 1030;
animal(ianim).chns = [4 40];

animal(ianim).tags(1).tag = 'CL';
animal(ianim).tags(1).expt_list = [101:105 107];

animal(ianim).tags(2).tag = 'OL';
animal(ianim).tags(2).expt_list = 201:202;

animal(ianim).tags(3).tag = 'gray';
animal(ianim).tags(3).expt_list = [106 109];

animal(ianim).tags(4).tag = 'passive';
animal(ianim).tags(4).expt_list = [101:105 107];

animal(ianim).tags(5).tag = 'active';
animal(ianim).tags(5).expt_list = [101:105 107];

ianim = ianim +1;
%%
animal(ianim).name = 'M130920_BALL';
animal(ianim).series = 1025;
animal(ianim).chns = [8 40];

animal(ianim).tags(1).tag = 'CL';
animal(ianim).tags(1).expt_list = 101:105;

animal(ianim).tags(2).tag = 'OL';
animal(ianim).tags(2).expt_list = 201:203;

animal(ianim).tags(3).tag = 'gray';
animal(ianim).tags(3).expt_list = [106 108];

animal(ianim).tags(4).tag = 'passive';
animal(ianim).tags(4).expt_list = 101:105;

animal(ianim).tags(5).tag = 'active';
animal(ianim).tags(5).expt_list = 101:105;

ianim = ianim +1;
%%
% animal(ianim).name = 'M131122_BALL';
% animal(ianim).series = 1212;
% animal(ianim).chns = [20 4];
% 
% animal(ianim).tags(1).tag = 'CL';
% animal(ianim).tags(1).expt_list = 107:110;
% 
% animal(ianim).tags(2).tag = 'OL';
% animal(ianim).tags(2).expt_list = 201;
% 
% animal(ianim).tags(3).tag = 'gray';
% animal(ianim).tags(3).expt_list = 107:110;

for ianim = 1:length(animal)
    for itag = 1:length(animal(ianim).tags)
        try
            es = VRLoadMultipleExpts(animal(ianim).name, animal(ianim).series, animal(ianim).tags(itag).expt_list,'SPIKES');
            save([animal(ianim).name '_' num2str(animal(ianim).series) '_' animal(ianim).tags(itag).tag], 'es');
%             powAnalysis(animal(ianim).name, animal(ianim).series, animal(ianim).tags(itag).expt_list, animal(ianim).tags(itag).tag, animal(ianim).chns, 1);
        catch
            display(['Could not run for animal ' animal(ianim).name ' and tag ' animal(ianim).tags(itag).tag]);
        end
    end
end

end