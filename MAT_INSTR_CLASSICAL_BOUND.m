rng('default');

% the last one is party D which is the channel
ins = [3,2,2,2];
outs = [2,2,2,2];

U = give_Joe_U();
BigU = kron(eye(2),U);
KrausU = {U};
ChoiU = ChoiMatrix(KrausU);


% % % iniP_proj = givePprojDET();%givePprojRANDgeneral(3,ins(1:end-1),outs(1:end-1));
% % % 
% % % U = give_Joe_U();
% % % BigU = kron(eye(2),U);
% % % KrausU = {U};
% % % ChoiU = ChoiMatrix(KrausU);
% % % 
for w = 1:ins(4)
    for d = 1:outs(4)
        %iniChannelsinst{w}{d} = ChoiU;%{give_Joe_U()};%{};
        iniChannelsinst{w}{d} = ChoiU;%{giveChannelRAND(2,4)};
    end
end

Pproj = givePprojRAND();
inichannel = iniChannelsinst;
output = criticalvisibilityInstrumental(Pproj, inichannel, ins, outs);
alpha = output{1};
bellcoeffs = output{2};
constr2 = output{4};
disp(join([string(alpha),'=',string(dispBellCoeffsINSTR(bellcoeffs,ins, outs))]));


out=SeeSawOverChannelInstr(alpha, bellcoeffs, Pproj, ins, outs);
disp(out{1})
disp(out{2})