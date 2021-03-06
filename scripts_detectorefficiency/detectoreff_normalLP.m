% add the functions folders to path
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);


ins = [3,2,2];
outs = [2,2,2];
nr_parties = length(ins);

load('best_for_sqrt3.mat');
p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs, ins, outs);
p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs, ins, outs);
[newalpha, newbellcoeffs, LPstatus, dual_alpha] = BroadcastLP(p_entangled, p_uniform);



% for eta=flip(0.545:0.05:0.9)
%    efficiencies = eta*[1,1,1];
%    p_dec_eff = giveDetectorEfficiencydistrib(p_entangled, efficiencies, ins, outs);
%    [bellcoeffs, LPstatus] = BroadcastLPfeasibility(p_dec_eff);
%    fprintf("detection efficiency: %f LPstatus= %d (1 is infeasible, 0 feasible)\n", eta, LPstatus);
% end

MAX_ITER_OUTER_LOOP = 10;
MAX_ITER_INNER_LOOP = 10;
ETA_CONVERGENCE_TOL = 1e-6;

list_of_conv_etas = [];
list_of_conv_channels = {};
list_of_conv_povms = {};

yalmip_optimizer = BroadcastLPfeasibility_optimizer(ins,outs+1);

%% Start with a good point
channel = {give_Joe_U()}; %best_channels{1};
POVMs = givePprojDET(); %best_povm{1};
p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs , ins, outs);
[eta, bellcoeffs] = BroadcastSlowLPdetectorefficiency(yalmip_optimizer,p_entangled,ins,outs,ETA_CONVERGENCE_TOL);


outer_iteration = 1;
while outer_iteration <= MAX_ITER_OUTER_LOOP
    fprintf("\nOuter loop iteration: %d\n", outer_iteration);

    inner_iteration = 1;
    change_eta = 100000;
    while inner_iteration < MAX_ITER_INNER_LOOP && change_eta > ETA_CONVERGENCE_TOL && abs(eta)>1e-2 && abs(eta-1)>1e-2
        old_eta = eta;
        [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(0), POVMs, channel);
        p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs, ins, outs);
        [new_eta, bellcoeffs] = BroadcastSlowLPdetectorefficiency(yalmip_optimizer,p_entangled,ins,outs,ETA_CONVERGENCE_TOL);
        change_eta = abs(new_eta-old_eta);
        fprintf("for inner_iteration=%d eta=%f",inner_iteration, new_eta);
        inner_iteration = inner_iteration + 1;
        if new_eta < old_eta 
            eta = new_eta;
            break;
        end
    end
    fprintf("\n eta is now: %f\n", eta);
    
    while abs(eta)<1e-2 && abs(eta-1)<1e-2
        POVMs = givePprojRANDgeneral(ins);
        channel = RandomSuperoperator([2,4]);
        p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs, ins, outs);
        [eta, bellcoeffs] = BroadcastSlowLPdetectorefficiency(yalmip_optimizer,p_entangled,ins,outs,ETA_CONVERGENCE_TOL);
        fprintf("Detector efficiency given the initial conditions: %f\n", eta);
    end
    
    
    list_of_conv_etas = [list_of_conv_etas, eta];
    list_of_conv_channels{outer_iteration} = channel;
    list_of_conv_povms{outer_iteration} = POVMs;
    outer_iteration = outer_iteration + 1;
end

ScenarioFilename = '322';
filename = strcat(ScenarioFilename,'best_detector_efficiency','.mat');
save(filename);
