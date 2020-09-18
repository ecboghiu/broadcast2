U = give_Joe_U();
BigU = kron(eye(2),U);
KrausU = {U};
ChoiU = ChoiMatrix(KrausU);

ins = {[1,2,3],[1,2],[1,2]};
outs = {[1,2],[1,2],[1,2]};


% out = MetaCompleteRound(ins, outs);
% 
% % Analyze the outputs
% alpha = out{1};
% channel = out{2};
% choi_form_channel = ChoiMatrix(channel);
% fprintf("final alpha: %f\n", alpha);
% bellcoeffs = out{4};
% Pproj = out{3};
% 
% %disp(choi_form_channel);
% %disp(ChoiU);
% 
% disp(join([string(alpha),'=',string(dispBellCoeffs(bellcoeffs,ins, outs))]));
% [expression, scaling] = dispBellCoeffsCorrelators(bellcoeffs,ins, outs);
% disp(join([string(vpa(alpha*scaling,3)),'=',string(expression)]));

a = LocalBoundOfIneq(bellcoeffs,ins,outs);
disp(a{1})
disp(a{2})

% Uncomment this for debugging
% state = final_state(ini_state(alpha),ChoiU);
% [newPproj,alpha] = SeeSawOverASingleParty(1, state, bellcoeffs, Pproj, ins, outs)
% disp(alpha);
% state = final_state(ini_state(0.8),{give_Joe_U()});
% Pproj = givePprojRAND();
% for party=1:3
%     for x=ins{party}
%         sull = 0;
%         for a=outs{party}
%             sull = sull + Pproj{party}{x}{a};
%         end
%         disp(sull);
%     end
% end

