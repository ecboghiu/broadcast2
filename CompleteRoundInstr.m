function out = CompleteRoundInstr(bellcoeffsFIRST, belllocalboundFIRST, ins, outs)

abschangealpha = 1e6;

alpha=1;
%while alpha>1-1e-2
    %iniP_proj = givePprojRAND2();
    iniP_proj = givePprojRANDgeneral(length(ins),ins,outs);
    %iniP_proj = givePprojRANDmaxEBI();
    %iniP_proj = givePprojDET();
    %iniChannel = RandomSuperoperator([2,4]);
    iniChannel = {giveChannelRAND(2,4)};
    %iniChannel = {give_Joe_U()};
    
    oldalpha=alpha;
    %outputcritvis = criticalvisibility(iniP_proj, iniChannel, ins, outs);
    
    %disp(iniChannel{1});
    %alpha = outputcritvis{1};
    %bellcoeffs = outputcritvis{2};
    fprintf("Alpha given initial conditions (should be 1): %f\n", alpha);
%end

ABS_TOL = 1e-6;
MAX_ITER = 1000;
iteration = 0;
abschange = 1e6;

Pproj = iniP_proj;
channel = iniChannel;
bellcoeffs = bellcoeffsFIRST;
belllocalbound = belllocalboundFIRST;
listofalphas = [];
while abschangealpha>ABS_TOL && iteration<MAX_ITER
    inistate = final_state(ini_state(alpha),channel);
    [newPproj,~,newchannel] = SeeSawOverAllParties(inistate, bellcoeffs, Pproj, ins, outs, alpha, channel);
    Pproj = newPproj;
    channel = newchannel;
    
    %disp(ChoiMatrix(channel));
    %[newchannel,newobjval2] = SeeSawOverChannel(alpha, bellcoeffs, newPproj, ins, outs);
    %channel = newchannel;
    
    %abschange = abs(newobjval2 - objval);
    %fprintf("alpha, bellPARTIES, bellPARTIES+CHANNEL: %f %f %f\n", alpha, newobjval1, newobjval2); 
    %objval = newobjval2;
    
    %DisplaySumPOVMS(Pproj,ins,outs);
    
    oldalpha = alpha;
    outputcritvis = criticalvisibilityFixedBell2(Pproj, channel, ins, outs);
    alpha1         = outputcritvis{1};
    %bellcoeffs    = outputcritvis{2};
    belllocalbound = ClassicalOptInequality2(bellcoeffs,ins,outs);
    disp(belllocalbound)
    
    alpha2 = VisibilityOfBellIneq(bellcoeffs, belllocalbound, Pproj, channel, ins, outs);
    
    listofalphas = [listofalphas, alpha1, alpha2];
    
    abschangealpha = abs(alpha-oldalpha);
    fprintf("\n alpha1 alpha2 =%f %f\n", alpha1, alpha2);
    
    iteration = iteration + 1;
end

out=cell(1);
out{1}=alpha;
out{2}=channel;
out{3}=Pproj;
out{4}=bellcoeffs;
out{5}=listofalphas;

end

