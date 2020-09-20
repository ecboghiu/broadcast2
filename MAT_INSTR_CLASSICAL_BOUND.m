rng('default');

% the last one is party D which is the channel
ins = [3,2,2,2];
outs = [2,2,2,2];

% % % iniP_proj = givePprojDET();%givePprojRANDgeneral(3,ins(1:end-1),outs(1:end-1));
% % % 
% % % U = give_Joe_U();
% % % BigU = kron(eye(2),U);
% % % KrausU = {U};
% % % ChoiU = ChoiMatrix(KrausU);
% % % 
iniChannelsinst = {{0}};
for w = 1:ins(4)
    for d = 1:outs(4)
        iniChannelsinst{w}{d} = ChoiU;%{give_Joe_U()};%{giveChannelRAND(2,4)};
        iniChannelsinst{w}{d} = ChoiU;
    end
end
% % % 
% % % 
% % % meas = iniP_proj;
% % % channels = iniChannels;
% % % out = criticalvisibilityINSTR(givePprojDET(),{{ChoiU}},ins,outs);
% % % bell1 = out{2};
% % % constr1 = out{5};
% % % disp(join([out{1},'=',string(dispBellCoeffsINSTR(bell1,ins,outs))]));

U = give_Joe_U();
BigU = kron(eye(2),U);
KrausU = {U};
ChoiU = ChoiMatrix(KrausU);

inichannel = iniChannelsinst;
output = criticalvisibility_old(givePprojDET(), inichannel, ins, outs);
alpha = output{1};
bellcoeffs = output{2};
constr2 = output{4};
disp(join([string(alpha),'=',string(dispBellCoeffsINSTR(bellcoeffs,ins, outs))]));