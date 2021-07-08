mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

%%
load('matlabworkspace_3221_2221.mat');

bellcoeffs = best_everything{:, 4};
mat = fromBellToCorrMat_INSTR(ins,outs);

new_bell = bellcoeffs(:)' * mat';
new_bell_corr = ToCorrelatorNotationINSTR_sym(new_bell, ins, outs);

[C,T] = coeffs(new_bell_corr);
C(abs(C)<1e-6) = 0;
min_coeff = min(abs(C));
C = C/min_coeff;
new_bell_corr2 = dot(C,T);

new_bell_corr2 = vpa(new_bell_corr2, 3);
disp(new_bell_corr2)