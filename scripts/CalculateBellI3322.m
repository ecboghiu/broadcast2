mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

% Fix the seed for debugging.
%rng(1234); % this seems to saturate the local bound
rng('default');

C1 = sym('C1');
C2 = sym('C2');
B1 = sym('B1');
B2 = sym('B2');
B3 = sym('B3');
B4 = sym('B4');
A1 = sym('A1');
A2 = sym('A2');
A3 = sym('A3');
A4 = sym('A4');
  
% WARNING for I3322 set FLAG_Use01obsInsteadOfCorrelator to true
% here Ai Bj mean observables in 0,1 outcomes
% see eq (3) in https://arxiv.org/pdf/1006.3032.pdf
ins = [3,3];
outs = [2,2];
I3322 = -A2 - B1 - 2*B2 + A1*B1 + A1*B2 + A2*B1 + A2*B2 ...
                        - A1*B3 + A2*B3 - A3*B1 + A3*B2;

I3322manualProj = -   A2*(sym('Pi_2_1_1')+sym('Pi_2_1_2')) ...
                  -   B1*(sym('Pi_1_2_1')+sym('Pi_1_2_2'))...
                  - 2*B2*(sym('Pi_1_1_1')+sym('Pi_1_1_2')) ...
                        + A1*B1 + A1*B2 + A2*B1 + A2*B2 ...
                        - A1*B3 + A2*B3 - A3*B1 + A3*B2;
I3322manualProj = subs(I3322manualProj, A1, sym('Pi_1_1_1'));        
I3322manualProj = subs(I3322manualProj, A2, sym('Pi_1_2_1'));  
I3322manualProj = subs(I3322manualProj, A3, sym('Pi_1_3_1'));  
I3322manualProj = subs(I3322manualProj, B1, sym('Pi_2_1_1'));  
I3322manualProj = subs(I3322manualProj, B2, sym('Pi_2_2_1'));  
I3322manualProj = subs(I3322manualProj, B3, sym('Pi_2_3_1'));  
I3322manualProj = expand(I3322manualProj);
I3322prob = FromProjectorNotationToProbability(I3322manualProj,ins,outs);
bellcoeffsI3322 = GetBellCoeffsFromProbSymIneq(I3322prob,ins,outs);
LBoundI3322 = ClassicalOptInequalityNoBroadcast(bellcoeffsI3322)
ins = [4,3,2];
outs = [2,2,2];
ModifiedI3322 = I3322manualProj * (C1 + C2) + LBoundI3322 * A4 * (C1 - C2);
ModifiedI3322 = subs(ModifiedI3322, C1, sym('Pi_3_1_1'));  
ModifiedI3322 = subs(ModifiedI3322, C2, sym('Pi_3_2_1'));  
ModifiedI3322 = subs(ModifiedI3322, A4, sym('Pi_1_4_1'));  
ModifiedI3322 = expand(ModifiedI3322);
ModifiedI3322prob = FromProjectorNotationToProbability(ModifiedI3322,ins,outs);
bellcoeffs = GetBellCoeffsFromProbSymIneq(ModifiedI3322prob,ins,outs);  
             
LBoundEBI = ClassicalOptInequalityNoBroadcast(bellcoeffs)

save('bellcoeffsI3322broadcast.mat','bellcoeffs','ins','outs');