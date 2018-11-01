function [arrayLayout, ArrayNo, shanks, tet] = MichiganGetLayout(animal, iseries, nchan, type)
% MichiganGetLayout - gives the physical arrangement of Michigan array channels in an NxN matrix
%
% [arrayLayout, shanks, tet] = MichiganGetLayout operates on
% current choice of animal and series (global PICK)
%
% [arrayLayout, shanks, tet] = MichiganGetLayout(animal, iseries) 
% lets you specify animal and series
%
% [arrayLayout, shanks, tet] = MichiganGetLayout(animal, iseries, type)
%
% It returns a matrix with the layout of the channels, and the number of
% shanks, whether or not in tetrode configuration and whether they are 
% tetrodes if tet = 1
%
% If the argument type is not provided, the function might
% prompt the user to specify the serial number of the array that was used in a
% particular series.
%
% 2012-01 MC turned it into a simple call to a method of object Probe.

global PICK

if nargin < 4
    type = [];
end

if nargin< 3
     nchan = [];
end

if nargin < 2
    animal = PICK.animal;
    iseries = PICK.iseries;
end

% warning('MichiganGetLayout is becoming obsolete...');

myProbe = Probe( animal, iseries, nchan, type );
arrayLayout = myProbe.Layout;
ArrayNo     = myProbe.Type;
shanks      = myProbe.nShanks;
tet         = ~isempty(myProbe.Tetrodes);




