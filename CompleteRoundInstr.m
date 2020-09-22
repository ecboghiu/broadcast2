function out = CompleteRoundInstr(bellcoeffsFIRST, belllocalboundFIRST, ins, outs)

abschangealpha = 1e6;

dimA=2;
dimB1=2;
dimB2=2;

alpha=1;
%while alpha>1-1e-2
    %iniP_proj = givePprojRAND2();
    iniP_proj = givePprojRANDgeneral(length(ins)-1,ins(1:end-1),outs(1:end-1));
    %iniP_proj = givePprojRANDmaxEBI();
    %iniP_proj = givePprojDET();
    %iniChannel = RandomSuperoperator([2,4]);

    U = give_Joe_U();
    %BigU = kron(eye(2),U);
    KrausU = {U};
    ChoiU = ChoiMatrix(KrausU);
    for w = 1:ins(4)
    %for d = 1:outs(4)
        %iniChannelsinst{w}{d} = ChoiU;%{give_Joe_U()};%{};
        iniChannels{w}= RandomPOVM(dimA*dimB1*dimB2,outs(4));%{ChoiU/2,ChoiU/2};%{giveChannelRAND(2,4)};
    %end
    end
    
    
    oldalpha=alpha;
    output = criticalvisibilityInstrumental(iniP_proj, iniChannels, ins, outs);
    alpha = output{1};
    inibellcoeffs = output{2};
    %constr2 = output{4};
    %disp(join([string(alpha),'=',string(dispBellCoeffsINSTR(bellcoeffs,ins, outs))]));
    
    %outputcritvis = criticalvisibility(iniP_proj, iniChannel, ins, outs);
    
    %disp(iniChannel{1});
    %alpha = outputcritvis{1};
    %bellcoeffs = outputcritvis{2};
    fprintf("Alpha given initial conditions (should be 1): %f\n", alpha);
%end

ABS_TOL = 1e-6;
MAX_ITER = 1000;
iteration = 1;
abschange = 1e6;

Pproj = iniP_proj;
channels = iniChannels;
bellcoeffs = inibellcoeffs;
belllocalbound = alpha;
listofalphas = {};
while abschangealpha>ABS_TOL && iteration<MAX_ITER
    %inistate = final_state(ini_state(alpha),channels);
    [newPproj,~,newchannels] = SeeSawOverAllPartiesInstrumental(bellcoeffs, Pproj, ins, outs, alpha, channels);
    Pproj = newPproj;
    channels = newchannels;
    
    %disp(ChoiMatrix(channel));
    %[newchannel,newobjval2] = SeeSawOverChannel(alpha, bellcoeffs, newPproj, ins, outs);
    %channel = newchannel;
    
    %abschange = abs(newobjval2 - objval);
    %fprintf("alpha, bellPARTIES, bellPARTIES+CHANNEL: %f %f %f\n", alpha, newobjval1, newobjval2); 
    %objval = newobjval2;
    
    %DisplaySumPOVMS(Pproj,ins,outs);
    
    oldalpha=alpha;
    output = criticalvisibilityInstrumental(iniP_proj, iniChannels, ins, outs);
    alpha = output{1};
    bellcoeffs = output{2};
    
    abschangealpha = abs(alpha-oldalpha);
    fprintf("alpha =%f %f\n", alpha);
    
    listofalphas{iteration} = {alpha, bellcoeffs, Pproj, channels};
    
    % TODO Fixed bell?
    %oldalpha = alpha;
    %outputcritvis = criticalvisibilityFixedBell2(Pproj, channels, ins, outs);
    %alpha         = outputcritvis{1};
    %bellcoeffs    = outputcritvis{2};
    
    %TODO be able to find local bound of bell ineq
    %belllocalbound = ClassicalOptInequality2(bellcoeffs,ins,outs);
    %disp(belllocalbound)
    
    %alpha2 = VisibilityOfBellIneq(bellcoeffs, belllocalbound, Pproj, channels, ins, outs);
    
    %listofalphas = [listofalphas, alpha1, alpha2];
    
    
    
    iteration = iteration + 1;
    
end

out=cell(1);
out{1}=alpha;
out{2}=channels;
out{3}=Pproj;
out{4}=bellcoeffs;
out{5}=listofalphas;

end

