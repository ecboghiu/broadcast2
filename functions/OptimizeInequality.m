function [outputArg1,outputArg2] = OptimizeInequality(bellcoeffs, ins, outs)
%nrparties = length(ins);

abschangealpha = 1e6;


%iniP_proj = givePprojRAND2();
iniP_proj = givePprojRAND();
%iniP_proj = givePprojDET();
%iniChannel = RandomSuperoperator([2,4]);
iniChannel = {giveChannelRAND(2,4)};
%iniChannel = {give_Joe_U()};


%disp(iniChannel{1});
alpha = 1;

ABS_TOL = 1e-6;
MAX_ITER = 1000;
iteration = 0;
abschange = 1e6;

Pproj = iniP_proj;
channel = iniChannel;
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
    outputcritvis = criticalvisibility(Pproj, channel, ins, outs);
    alpha         = outputcritvis{1};
    bellcoeffs    = outputcritvis{2};
    
    abschangealpha = abs(alpha-oldalpha);
    fprintf("\n alpha=%f, abschangealpha=%g\n\n", alpha, abschangealpha);
    
    iteration = iteration + 1;
end

out=cell(1);
out{1}=alpha;
out{2}=channel;
out{3}=Pproj;
out{4}=bellcoeffs;

end

