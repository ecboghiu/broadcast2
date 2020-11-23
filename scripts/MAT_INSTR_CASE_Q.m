rng('default');

% the last one is party D which is the channel
ins = [3,2,2,1];
outs = [2,2,2,1];

out = MetaCompleteRoundINSTR(ins, outs);
listalphas = out{5};
alpha=out{1};
channels=out{2};
Pproj=out{3};
bellcoeffs=out{4};
listofalphas=out{5};
