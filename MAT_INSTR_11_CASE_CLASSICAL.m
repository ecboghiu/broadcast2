rng('default');

% the last one is party D which is the channel
ins = [3,2,2,1];
outs = [2,2,2,1];

U = give_Joe_U();
BigU = kron(eye(2),U);
KrausU = {U};
ChoiU = ChoiMatrix(KrausU);

dimA = 2;
dimB1 = 2;
dimB2 = 2;
for w = 1:ins(4)
    %for d = 1:outs(4)
        iniChannelsinst{w} = {ChoiU/2,ChoiU/2};%{give_Joe_U()};%{};
        %iniChannelsinst{w} = RandomPOVM(dimA*dimB1*dimB2,outs(4));
    %end
end

Pproj = givePprojDET();
inichannels = iniChannelsinst;
output = criticalvisibilityInstrumental(Pproj, inichannels, ins, outs);
alpha = output{1};
bellcoeffs = output{2}; bellcoeffs(abs(bellcoeffs)<1e-6)=0;
constr2 = output{4};
disp(join([string(alpha),'=',string(dispBellCoeffsINSTR(bellcoeffs,ins, outs))]));