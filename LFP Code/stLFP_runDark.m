%%
animal = 'M130412_BALL';
iseries = 514;
iexp  = 105;
iexp2 = 106;

% [a412_e1.es a412_e1.spikeTimes a412_e1.lfp] = stLFP_oldData(animal, iseries, iexp, 4, 1);
[a412_e2.es a412_e2.spikeTimes a412_e2.lfp] = stLFP_oldData(animal, iseries, iexp2, 4, 1);

% save('a412_e1', 'a412_e1', '-v7.3')
save('a412_e2', 'a412_e2', '-v7.3')
clear

%%
animal = 'M130503_BALL';
iseries = 520;
iexp  = 106;
iexp2 = 107;

[a503_e1.es a503_e1.spikeTimes a503_e1.lfp] = stLFP_oldData(animal, iseries, iexp, 4, 1);
[a503_e2.es a503_e2.spikeTimes a503_e2.lfp] = stLFP_oldData(animal, iseries, iexp2, 4, 1);

save('a503_e1', 'a503_e1', '-v7.3')
save('a503_e2', 'a503_e2', '-v7.3')
clear

%%
animal = 'M130504_BALL';
iseries = 524;
iexp  = 105;
iexp2 = 106;

[a504_e2.es a504_e2.spikeTimes a504_e2.lfp] = stLFP_oldData(animal, iseries, iexp2, 4, 1);
[a504_e1.es a504_e1.spikeTimes a504_e1.lfp] = stLFP_oldData(animal, iseries, iexp, 4, 1);

save('a504_e1', 'a504_e1', '-v7.3')
save('a504_e2', 'a504_e2', '-v7.3')
clear