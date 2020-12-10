%% Working directories
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

%% Fix the scenario
party_ins = [3,2,2];
party_outs = [2,2,2];
instr_ins = [2];
instr_outs = [2];
dims_in = [2];
dims_out = [4];
ins = [party_ins, instr_ins];
outs = [party_outs, instr_outs];

%% To print the scenario.
aux_ins = string([party_ins, instr_ins]);
aux_outs = string([party_outs, instr_outs]);
Scenario=strcat(aux_ins{:},'-',aux_outs{:});
fprintf("Instrumental scenario = xyzw-abcd (w instrument input, d instrument output) = %s\n", Scenario);
diaryname=strcat('mydiary',Scenario,'.txt');
diary(diaryname);

%% Fix the rng seed
semilla = sum(100*clock);
%semilla = 216391; % Fix it for debugging
%semilla = 1;
fprintf("Fixing random seed = %g\n", uint32(semilla));
rng(semilla,'twister');

%% How many total rounds
MAX_ITER_META = 1000;


%% The MAIN loop
final_round_alpha = [];
final_round_povm = {};
final_round_channels = {};
localbound = nan;

latest_alpha_meta = 0;
best_alpha = 0;
best_everything = {};
best_index = 1;
meta_iteration = 1;
while meta_iteration < MAX_ITER_META
    fprintf("\nRound %d.\n", meta_iteration);
    % While the visibility is close to 1 or above 1, keep changing the 
    % initial conditions until we get a behaviour outside the broadcast-local
    % polytope.
    alpha = 0;
    ALHPA_TOL_DIST_TO_1 = 1e-3;
    LPstatus = 0;
    %% Loop through initial conditions until they are good
    while abs(alpha-0)<ALHPA_TOL_DIST_TO_1 || abs(alpha-1)<ALHPA_TOL_DIST_TO_1 || alpha < 0.2 || LPstatus~=0 
        %% Initialize the POVMs.
        %iniPovms = givePprojRAND2();  % This uses QETLAB's RandomPOVM()
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

        if ~checkPOVMsAreGood(iniPovms,party_ins,party_outs)
           warning('Problem with POVMs\n');
        end
        if ~checkThatInstrChannelsAreGood(iniChannel, instr_ins, instr_outs, dims_in, dims_out)
           warning('Problem with the channel\n');
        end

        
        %% Run the LP with broadcast-local constraints.
        p_entangled = ProbMultidimArrayInstrumental(NoisyWernerState(0), iniPovms, iniChannel, ins, outs);
        p_uniform   = ProbMultidimArrayInstrumental(NoisyWernerState(1), iniPovms, iniChannel, ins, outs);
        [alpha, bellcoeffs, LPstatus] = BroadcastInstrumentLP(p_entangled, p_uniform, ins, outs);
        
        if LPstatus ~= 0
            cleaning_tol = 1e-12;
            iter = 1;
            while LPstatus ~= 0 && iter < 12
               fprintf("Using clean(p,%g). Previous LPstatus: %d\n", cleaning_tol, LPStatus);
               [final_alpha,bellcoeffs,LPstatus] = BroadcastInstrumentLP( ...
                   clean(p_entangled,cleaning_tol), ...
                   clean(p_uniform,cleaning_tol), ...
                   ins, outs);
               cleaning_tol = cleaning_tol * 10;
               iter = iter + 1;
            end
        end
        
        
        %[alpha, bellcoeffs] = BroadcastSlowLP(p_entangled, p_uniform, ins, outs);
        fprintf("Visibility given initial measurements and channels: %f\n", alpha);
        %localbound = ClassicalOptInequality_fromLPBroadcast_INSTR(bellcoeffs, ins, outs);

%         aux_bpent = bellcoeffs .* (p_entangled);
%         aux_bpuni = bellcoeffs .* (p_uniform);
%         aux_bdiff = bellcoeffs .* (p_entangled-p_uniform);
%         fprintf("With ini meas/channel: s�p1=%f, s�p2=%f, s�(p1-p2)=%f, localbound=%f\n", ...
%             sum(aux_bpent(:)), ...
%             sum(aux_bpuni(:)), ...
%             sum(aux_bdiff(:)), ...
%             localbound);
        
        %fprintf("Function check, alpha_vis_localbound=%f alpha_LP=%f\n", ...
        %    visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform), ...
        %    alpha);
    end

    %% Given initial conditions with ineq, do SDP and then LP again etc.
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

        [POVMs,finalObj,channels] = SeeSawOverAllPartiesInstrumental(bellcoeffs, NoisyWernerState(0), POVMs, channels, ins, outs);
        if ~checkPOVMsAreGood(POVMs,party_ins,party_outs)
           warning('Problem with POVMs\n');
        end
        if ~checkThatInstrChannelsAreGood(channels, instr_ins, instr_outs, dims_in, dims_out)
           warning('Problem with the channel\n');
        end

        oldalpha = alpha;
        p_entangled = ProbMultidimArrayInstrumental(NoisyWernerState(0), POVMs, channels, ins, outs);
        p_uniform   = ProbMultidimArrayInstrumental(NoisyWernerState(1), POVMs, channels, ins, outs);
 
        flag_some_prob_not_normalized = false;
        tol_prob_normalization = 1e-7;
        aux_ins_coords = ind2subv(ins, 1:prod(ins));
        for aux=1:size(aux_ins_coords,1)
           aux_ins_cell=num2cell(aux_ins_coords(aux,:));
           probvec = p_entangled(aux_ins_cell{:},:,:,:,:);
           if abs(sum(probvec(:))-1)>tol_prob_normalization
               flag_some_prob_not_normalized = true;
               pvec2 = p_uniform(aux_ins_cell{:},:,:,:,:);
               warning("Probability not normalized to precision. Example sum(prob): %g %g\n,", sum(probvec(:)), sum(pvec2(:)));
               break;
           end
        end
        
        [newalpha, bellcoeffs, LPstatus] = BroadcastInstrumentLP(p_entangled, p_uniform, ins, outs);
        if LPstatus ~= 0
            cleaning_tol = 1e-12;
            iter = 1;
            while LPstatus ~= 0 && iter < 9 && newalpha > 0.7
               fprintf("Using clean(p,%g). Previous LPstatus: %d\n", cleaning_tol, LPstatus);
               [newalpha,bellcoeffs,LPstatus] = BroadcastInstrumentLP( ...
                   clean(p_entangled,cleaning_tol), ...
                   clean(p_uniform,cleaning_tol), ...
                   ins, outs);
               cleaning_tol = cleaning_tol * 10;
               iter = iter + 1;
            end
            fprintf("If cleaning didn't help, check the validity of the POVMs and channels.\n");
            fprintf("Now applying bisection algorithm.\n");
            [probablybad_newalpha,bellcoeffs] = BroadcastSlowLP( ...
                   p_entangled, ...
                   p_uniform, ...
                   ins, outs);
            if 0 < probablybad_newalpha && probablybad_newalpha < 0.8
               LPstatus = 0; 
            end
            newalpha = 0; %% TODO FIX this, sometimes the alpha I get is bad...
        end
        
%         aux_bpent = bellcoeffs .* (p_entangled);
%         aux_bpuni = bellcoeffs .* (p_uniform);
%         aux_bdiff = bellcoeffs .* (p_entangled-p_uniform);
%         fprintf("With optimized meas/channel: s�p1=%f, s�p2=%f, s�(p1-p2)=%f, localbound=%f, alpha=%f, LPstatus=%f\n", ...
%             sum(aux_bpent(:)), ...
%             sum(aux_bpuni(:)), ...
%             sum(aux_bdiff(:)), ...
%             localbound, ...
%             newalpha, ...
%             LPstatus);

        %fprintf("Function check, alpha1=%f alpha2=%f diff=%f\n", ...
        %    visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform), ...
        %    newalpha,  ...
        %    visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform) - newalpha);
        
        %fprintf("Visibility after optimizing: %f\n", visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform));
        if LPstatus ~= 0 || newalpha >= 1-1e-3
            %disp('LP not solved correctly');
            warning("LP not solved correctly. Trying another set of initial points. LpStatus=%f, alpha=%f\n", LPstatus, newalpha);
            break;
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
    end
    if ~isempty(list_of_alphas)
        if alpha > best_alpha
           best_alpha = alpha;
           best_povm =  list_of_povms{end};
           best_channels = list_of_channels{end};
           best_bellcoeffs = list_of_bellcoeffs{end};
           best_everything{best_index} = {best_alpha, best_povm, best_channels, best_bellcoeffs};
           best_index = best_index + 1;
        end
        fprintf("Best visibility: %f\n", best_alpha);
    end
    
    meta_iteration = meta_iteration + 1;
end

save(strcat('matlabworkspace_asini',Scenario,'.mat'));
