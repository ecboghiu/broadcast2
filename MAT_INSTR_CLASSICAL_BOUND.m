rng('default');

% the last one is party D which is the channel
ins = [3,2,2,1];
outs = [2,2,2,1];

% U = give_Joe_U();
% BigU = kron(eye(2),U);
% KrausU = {U};
% ChoiU = ChoiMatrix(KrausU);


% % % iniP_proj = givePprojDET();%givePprojRANDgeneral(3,ins(1:end-1),outs(1:end-1));
% % % 
% % % U = give_Joe_U();
% % % BigU = kron(eye(2),U);
% % % KrausU = {U};
% % % ChoiU = ChoiMatrix(KrausU);
% % % 

% % for w = 1:ins(4)
% %     %for d = 1:outs(4)
% %         %iniChannelsinst{w}{d} = ChoiU;%{give_Joe_U()};%{};
% %         iniChannelsinst{w} = {ChoiU/2,ChoiU/2};
% %     %end
% % end
% % 
% % Pproj = givePprojRAND();
% % inichannels = iniChannelsinst;
% % output = criticalvisibilityInstrumental(Pproj, inichannels, ins, outs);
% % alpha = output{1};
% % bellcoeffs = output{2};
% % constr2 = output{4};
% % disp(join([string(alpha),'=',string(dispBellCoeffsINSTR(bellcoeffs,ins, outs))]));
% % 
% % 
% % out=SeeSawOverAllPartiesInstrumental(bellcoeffs, Pproj, ins, outs, alpha, inichannels);
% % disp(out{1})
% % disp(out{2})

 out = MetaCompleteRoundINSTR(ins, outs);
 listalphas = out{5};
alpha=out{1};
channels=out{2};
Pproj=out{3};
bellcoeffs=out{4};
listofalphas=out{5};
