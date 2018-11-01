function protocol = ProtocolInitialize(xfilename,quiet)
% ProtocolInitialize initializes a protocol data structure
%
% protocol = ProtocolInitialize
%
% protocol = ProtocolInitialize(xfilename)
%
% protocol = ProtocolInitialize(xfilename,'quiet') suppresses any text
% output (important e.g. when called from a Visual Basic program).

% part of Spikes
%
% 2000-09 Matteo Carandini
% 2006-08 MC updated and extended so one can declare an xfile

% To see what is in a Protocol, do load('Z:\Data\trodes\CATZ009\19\1\Protocol');

if nargin<2
    quiet = 'loud';
end

if nargin< 1
    xfilename = [];
end

    protocol.xfile      = '';
    protocol.adapt.flag = 0;
    protocol.nstim      = [];
    protocol.npfilestimuli = [];
    protocol.npars      = [];
    protocol.pars       = [];
    protocol.parnames   = {};
    protocol.pardefs 	= {};
 
    protocol.animal     = '';
    protocol.iseries	= [];
    protocol.iexp       = [];
    protocol.nrepeats   = 0;
    protocol.seqnums    = 0;
          
    % These are not key attributes of the Protocol file:
%     protocol.blankstims 	= [];
%     protocol.blankpars      = [];
%     protocol.activepars 	= [];
%     protocol.description 	= [];

if isempty(xfilename)
    return
end

%  read the xfile

x = XFileLoad( xfilename, quiet );

protocol.xfile      = x.name;
protocol.parnames   = x.parnames;
protocol.pardefs    = x.pardefs;
protocol.npars      = x.npars;



