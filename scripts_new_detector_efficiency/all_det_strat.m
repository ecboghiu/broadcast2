clear all

% add the functions folders to path
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

%% Exact ineq from paper
% TODO: Make this work!!!
% C1 = sym('C1');
% C2 = sym('C2');
% B1 = sym('B1');
% B2 = sym('B2');
% B3 = sym('B3');
% B4 = sym('B4');
% A1 = sym('A1');
% A2 = sym('A2');
% A3 = sym('A3');
% A4 = sym('A4');
% 
% ins = [3,2,2];
% outs = [2,2,2];
% PAPER = A1*B1*C1 + A1*B2*C2 + A2*B2*C2 - A2*B1*C1 + A1*B1*C2 + A1*B2*C1 + A2*B1*C2 - A2*B2*C1 - 2*A3*B1 + 2*A3*B2;
% 
% FLAG_Use01obsInsteadOfCorrelator = false;
% PaperProjector = ToProjectorNotation(PAPER,ins,outs,FLAG_Use01obsInsteadOfCorrelator);
% PaperProb = FromProjectorNotationToProbability(PaperProjector, ins, outs);
% PaperBellcoeffs = GetBellCoeffsFromProbSymIneq(PaperProb,ins,outs);
% PaperLocalbound =  ClassicalOptInequalityNoBroadcast(PaperBellcoeffs);
% bellcoeffs = PaperBellcoeffs;
% channel = {give_Joe_U()};
% POVMs = givePprojDET();


%%
load('goodJoeNumerical.mat');  % loads 'bellcoeffs' 'channel' 'POVMs'

ins = [3, 2, 2];
outs = [2, 2, 2];

p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs, ins, outs);
p_uniform = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs, ins, outs);
sum(bellcoeffs.*p_uniform,'all')

