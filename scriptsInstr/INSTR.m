mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

party_ins = [3,2,2];
party_outs = [2,2,2];
instr_ins = [1];
instr_outs = [1];
dims_in = [2];
dims_out = [4];
ins = [party_ins, instr_ins];
outs = [party_outs, instr_outs];


aux_ins = string([party_ins, instr_ins]);
aux_outs = string([party_outs, instr_outs]);
Scenario=strcat(aux_ins{:},'-',aux_outs{:});
%Scenario = '3223-2222';
fprintf("Instrumental scenario = xyzw-abcd (w instrument input, d instrument output) = %s\n", Scenario);
diaryname=strcat('mydiary',Scenario,'.txt');
diary(diaryname);

semilla = sum(100*clock);
semilla = 208887;
fprintf("Fixing random seed = %d\n", uint32(semilla));
rng(semilla,'twister');

MAX_ITER_META = 1000;

final_round_alpha = [];
final_round_povm = {};
final_round_channels = {};

latest_alpha_meta = 0;

best_alpha = 0;

meta_iteration = 1;
while meta_iteration < MAX_ITER_META
    fprintf("\nRound %d.\n", meta_iteration);
    % While the visibility is close to 1 or above 1, keep changing the 
    % initial conditions until we get a behaviour outside the broadcast-local
    % polytope.
    alpha = 0;
    ALHPA_TOL_DIST_TO_1 = 1e-3;
    LPstatus = 0;
    while abs(alpha-0)<ALHPA_TOL_DIST_TO_1 || LPstatus~=0
        %% Initialize the POVMs.
        %iniP_proj = givePprojRAND2();  % This uses QETLAB's RandomPOVM()
        iniPovms = givePprojRANDgeneral(party_ins);  % This generates random directions -as 
                                     % many as inputs- on the Bloch sphere
                                     % and gives the projectors onto these 
                                     % vectors. Assuming two outputs for
                                     % all parties.


        %% Choose an initial chanenl B -> B_1 $\otimes$ B_2.
        %iniChannel = RandomSuperoperator([2,4]);  % Using QETLAB's random
                                                   % superoperator
        %iniChannel = giveINSTRRandomSuperoperatorFamily(instr_ins, instr_outs, dims_in, dims_out);
        iniChannel = giveINSTRRandomIsometryFamily(instr_ins, instr_outs, dims_in, dims_out);
        %iniChannel = {give_Joe_U()};  % The optimal choice from the paper

        assert(checkPOVMsAreGood(iniPovms,party_ins,party_outs), 'Problem with POVMs');
        assert(checkThatInstrChannelsAreGood(iniChannel, instr_ins, instr_outs, dims_in, dims_out), 'Problem with the channel');
        
        %% Run the LP with broadcast-local constraints.
        p_entangled = ProbMultidimArrayInstrumental(NoisyWernerState(0), iniPovms, iniChannel);
        p_uniform   = ProbMultidimArrayInstrumental(NoisyWernerState(1), iniPovms, iniChannel);
        [alpha, bellcoeffs, LPstatus, dual_alpha] = BroadcastInstrumentLP(p_entangled, p_uniform, ins, outs);
        fprintf("Visibility given initial measurements and channels: %f\n", alpha);
        
        localbound = ClassicalOptInequality_fromLPBroadcast_INSTR(bellcoeffs, ins, outs);
        aux_bpent = bellcoeffs .* (p_entangled);
        aux_bpuni = bellcoeffs .* (p_uniform);
        aux_bdiff = bellcoeffs .* (p_entangled-p_uniform);
        fprintf("With ini meas/channel: s·p1, s·p2, s·(p1-p2), localbound: %f, %f, %f, %f\n", ...
            sum(aux_bpent(:)), ...
            sum(aux_bpuni(:)), ...
            sum(aux_bdiff(:)), ...
            localbound);

        %fprintf("Visibility after optimizing: %f\n", visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform));
    end

    % Parameters for stopping conditions for the next loop.
    ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
    MAX_ITER = 500;
    AbsDeltaAlpha = 1e6;
    iteration = 1;
    alpha_stdev = 1;

    POVMs = iniPovms;
    channels = iniChannel;
    list_of_alphas = [];
    list_of_povms = {};
    list_of_channels = {};
    list_of_bellcoeffs = {};
    while (AbsDeltaAlpha>ALPHA_CONVERGENCE_THRESHOLD && ...
            iteration<MAX_ITER && ...
            abs(alpha-0)>ALHPA_TOL_DIST_TO_1 && ...
            alpha_stdev > 1e-3)
        %[localbound, LPstatus] = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
        %fprintf("Local bound of bell inequality: %f\n", localbound);

        %[POVMs,finalObj,channels] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(alpha), POVMs, channels);
        [POVMs,finalObj,channels] = SeeSawOverAllPartiesInstrumental(bellcoeffs, NoisyWernerState(0), POVMs, channels, ins, outs);
        assert(checkPOVMsAreGood(POVMs,party_ins,party_outs),'Problem with POVMs');
        assert(checkThatInstrChannelsAreGood(channels, instr_ins, instr_outs, dims_in, dims_out), 'Problem with the channel'); % TODO FIX {instr}
        

        %outputcritvis = criticalvisibility_std(Pproj, channel, ins, outs);
        oldalpha = alpha;
        p_entangled = ProbMultidimArrayInstrumental(NoisyWernerState(0), POVMs, channels);
        p_uniform   = ProbMultidimArrayInstrumental(NoisyWernerState(1), POVMs, channels);
 
        aux_ins_coords = ind2subv(ins, 1:prod(ins));
        for aux=1:size(aux_ins_coords,1)
           aux_ins_cell=num2cell(aux_ins_coords(aux,:));
           probvec = p_entangled(aux_ins_cell{:},:,:,:,:);
           assert(abs(sum(probvec(:))-1)<1e-5,"Probability not normalized");
        end
        
        [newalpha, bellcoeffs, LPstatus, dual_alpha] = BroadcastInstrumentLP(p_entangled, p_uniform, ins, outs);
        localbound = ClassicalOptInequality_fromLPBroadcast_INSTR(bellcoeffs, ins, outs);
        if newalpha > 0.7
           error("Check wtf is wrong here"); 
        end
        
        aux_bpent = bellcoeffs .* (p_entangled);
        aux_bpuni = bellcoeffs .* (p_uniform);
        aux_bdiff = bellcoeffs .* (p_entangled-p_uniform);
        fprintf("With optimized meas/channel: s·p1, s·p2, s·(p1-p2), localbound, alpha, LPstatus: %f, %f, %f, %f, %f, %d\n", ...
            sum(aux_bpent(:)), ...
            sum(aux_bpuni(:)), ...
            sum(aux_bdiff(:)), ...
            localbound, ...
            newalpha, ...
            LPstatus);
        %fprintf("Visibility after optimizing: %f\n", visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform));
        if LPstatus ~= 0 || newalpha >= 1-1e-3
            %disp('LP not solved correctly');
            fprintf("LP not solved correctly. Trying another set of initial points.\n");
            break;
            %error('Check why infeasible problem.'); 
            list_of_alphas = [list_of_alphas, 0];
            list_of_povms{iteration} = 0;
            list_of_channels{iteration} = 0;
            list_of_bellcoeffs{iteration} = 0;
        else
            alpha = newalpha;
            list_of_alphas = [list_of_alphas, alpha];
            if max(size(list_of_alphas)) == 2 
                alpha_stdev = 1;
            elseif max(size(list_of_alphas)) == 1
                alpha_stdev = 1;
            else
                alpha_stdev = std(list_of_alphas(end-2:end));
            end
            list_of_povms{iteration} = POVMs;
            list_of_channels{iteration} = channels;
            list_of_bellcoeffs{iteration} = bellcoeffs;
            iteration = iteration + 1;
                
            AbsDeltaAlpha = abs(alpha-oldalpha);

            fprintf("Visibility after optimizing over POVMs and channels: %f\n", alpha);
        end
        %alpha = outputcritvis{1} + 0.01;
    end
    if ~isempty(list_of_alphas)
        if alpha > best_alpha
           best_alpha = alpha;
           best_povm =  list_of_povms{end};
           best_channels = list_of_channels{end};
           best_bellcoeffs = list_of_bellcoeffs{end};
        end
        %final_round_alpha = [final_round_alpha, list_of_alphas(end)];
        %final_round_povm{length(final_round_alpha)} = list_of_povms{end};
        %final_round_channels{length(final_round_alpha)} = list_of_channels{end};
        %latest_alpha_meta = final_round_alpha(end);
        fprintf("Best visibility: %f\n", best_alpha);
    end
    meta_iteration = meta_iteration + 1;
    
end

%best_alpha = max(final_round_alpha);
%index = find(final_round_alpha == best_alpha);
%best_povm = final_round_povm(index);
%best_channels = final_round_channels(index);

save(strcat('matlabworkspace_asini',Scenario,'.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%save 'data_optimal0.577Standard'
% 
% LBoundOurEBI = ClassicalOptInequality2(bellcoeffs,ins,outs);
% fprintf('LBoundOurEBI: %f\n', LBoundOurEBI);
% 
% inputstate = final_state(ini_state(alpha), channel);
% mixingstate = eye(size(inputstate))/size(inputstate,1);
% val=evaluate_bell_ineq_INSTR(bellcoeffs, 0, inputstate, Pproj, ins, outs);
% val
% 
% vis=VisibilityOfBellIneqWithNoChannel(bellcoeffs, LBoundOurEBI, Pproj, inputstate, mixingstate, ins, outs);
% vis
% 
% correlatorineq = dispBellCoeffsCorrelators(bellcoeffs,ins,outs);
% correlatorineq = simplify(vpa(correlatorineq,6));
% 
% probineq = ToProbabilityNotationIneqSym(correlatorineq,ins,outs);
% [bellcoeffsagain, ineqconstant] = GetBellCoeffsFromProbSymIneq(probineq,ins,outs);
% 
% [C,T] = coeffs(correlatorineq);
% 
% ineqconstant+evaluate_bell_ineq_INSTR(bellcoeffsagain, ineqconstant, final_state(ini_state(0.577), channel), Pproj, ins, outs)
% 
% LBoundOurEBI2 = ClassicalOptInequality2(bellcoeffsagain,ins,outs)
% fprintf("local bounds before and after: %f %f\n", LBoundOurEBI, LBoundOurEBI2);
