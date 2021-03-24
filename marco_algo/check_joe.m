diary arxiv1112_2626.txt

clear all

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
addpath(newdir);
addpath(newdir2);

rng('shuffle');

ins = [3,2,2];
outs = [2,2,2];

load('Joe0577.mat'); % load JoeCh JoeBell JoePovms

%%

bellcoeffs = JoeBell;
localbound_broadcast = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
channel = JoeCh;
POVMs = JoePovm;

p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs);
p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs);
alpha = visibilityOfBellInequality(bellcoeffs, localbound_broadcast, p_entangled, p_uniform);

%% Loop parameters
ALPHA_INI_ITER = 1000; % How many initial conditions to try for the initial point
MAX_ITER_INNER_LOOP = 200; % How many times to try the whole process
ALHPA_TOL_DIST_TO_POINT = 1e-3; % Util for excluding "bad" points such as visibility=0 or visibility=1 (if its zero, then most likely its an infeasibility)
ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
MAX_ITER_INNER = 50;
MAX_ITER_OUTER = 50;
INITIAL_VISIBILITY = 0;  % !! I'm using the convention of (1-p)psi+p*id instead of p*psi+(1-p)*id
CONST_CHI = 0.1;
DELTA_STATE_VIS = 0.001;

% Store the values after the program stops at a local optimum for all
% runs
inner_list_of_alphas = zeros(1,MAX_ITER_OUTER);
inner_list_of_povms = cell(1,MAX_ITER_OUTER);
inner_list_of_channels = cell(1,MAX_ITER_OUTER);

%outer_list_of_alphas = [];
%outer_list_of_povms = {};
%outer_list_of_channels = {};

% Store the best of teh local optimums
final_round_alpha = [];
final_round_povm = {};
final_round_channels = {};

NR_OF_INEQS=1;
ineq_nr=1;

iteration_inner_loop = 1;
while iteration_inner_loop <= MAX_ITER_INNER_LOOP 
    alpha=0;
    alpha_iteration=0;
    LPstatus = 0;

    visibility = INITIAL_VISIBILITY;
    %% FIND A GOOD INITIAL POINT
    while (abs(alpha-0)<ALHPA_TOL_DIST_TO_POINT || abs(alpha-1)<ALHPA_TOL_DIST_TO_POINT) && (alpha_iteration < ALPHA_INI_ITER)
        %POVMs = givePprojRANDgeneral(ins);  % Random projective measurements
        POVMs = givePprojRAND2(ins, outs); % QETLAB's random superopertor
        %channel = {giveChannelRAND(2,4)};  % Random isometry
        channel = RandomSuperoperator([2,4]);
        [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(visibility), POVMs, channel);
        p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs);
        p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs);
        alpha = visibilityOfBellInequality(bellcoeffs, localbound_broadcast, p_entangled, p_uniform);

        aux1 = bellcoeffs.*p_entangled; aux2 = bellcoeffs.*p_uniform;
        fprintf("Vis of ineq given random ini = %f. vis_state=%f Bell·p_ent=%f Bell·p_unif=%f bound_l=%f Ineq.nr.:%d/%d Ini.iter.:%d/%d\n", ...
            alpha, visibility, sum(aux1(:)), sum(aux2(:)), localbound_broadcast, ineq_nr, NR_OF_INEQS, alpha_iteration, ALPHA_INI_ITER);
        if abs(alpha-0)> ALHPA_TOL_DIST_TO_POINT && abs(alpha-1)> ALHPA_TOL_DIST_TO_POINT && alpha < 0.1
            fprintf("\nFound good initial condition. Bell value:%f NS2Bound: %f alpha=%f state_vis=%f\n", finalObj, localbound_broadcast, alpha, visibility);
            break;
        end
        alpha_iteration = alpha_iteration + 1;
    end
    fprintf("Given this point, optimize the visibility of the inequality.\n\n");
    %% DO SEE SAW STARTING FROM THIS POINT
    iteration_outer_loop = 1;
    old_visibility = 0;
    AbsDeltaStateVis = 1e6;
    while iteration_outer_loop <= MAX_ITER_OUTER && AbsDeltaStateVis > 1e-4
        AbsDeltaAlpha = 1e6;
        iteration_inner_loop = 1;
        % First optimize iteratively until the visibility converges to
        % something
        while (AbsDeltaAlpha > ALPHA_CONVERGENCE_THRESHOLD && ...
                iteration_inner_loop < MAX_ITER_INNER && ...
                abs(alpha-0) > ALHPA_TOL_DIST_TO_POINT)
            oldalpha = alpha;
            [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(visibility), POVMs, channel);
            p_entangled = ProbMultidimArray(final_state(NoisyWernerState(visibility), channel), POVMs);
            p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs);
            alpha = visibilityOfBellInequality(bellcoeffs, localbound_broadcast, p_entangled, p_uniform);
            AbsDeltaAlpha = abs(alpha-oldalpha);
            aux1 = bellcoeffs.*p_entangled; aux2 = bellcoeffs.*p_uniform;
            fprintf("Visibility of ineq %d/%d after optimizing = %f. vis_state=%f  Bellvalue=%f Bell·p_ent=%f Bell·p_unif=%f bound_l=%f\n", ...
            ineq_nr, NR_OF_INEQS, alpha, visibility, finalObj, sum(aux1(:)), sum(aux2(:)), localbound_broadcast);
        end

        inner_list_of_alphas = [inner_list_of_alphas, alpha];
        inner_list_of_povms{iteration_outer_loop} = POVMs;
        inner_list_of_channels{iteration_outer_loop} = channel;
%             
%             if alpha-(1-0.683)>1e-8 && abs(alpha-1)>1e-3 && abs(alpha)>1e-3
%                fprintf("\n\n\nFound what we're looking for. Bell value:%f NS2 Bound: %f alpha=%f\n\n\n\n", finalObj, localboundNS2,alpha);
%                save('final_workspace.mat');
%                error("Finishing program.");
%             end
        if abs(alpha)<ALHPA_TOL_DIST_TO_POINT
           fprintf("\nIneq vis converged to 0, reached local optimum. Restart with different initial conditions.\n\n");
           break;
        end
        % Then we go back a bit and decrease the noise and try again
        old_visibility = visibility;
        fprintf("\nIneq vis converged to: %f. Update state vis from %f to %f (%f below the local minima). DeltaStateVis=%f\n", alpha, visibility, visibility + alpha - DELTA_STATE_VIS, DELTA_STATE_VIS, + alpha - DELTA_STATE_VIS);
        visibility = visibility + alpha - DELTA_STATE_VIS; % remember I use the opposite convention for visibility
        AbsDeltaStateVis = abs(visibility-old_visibility);
        iteration_outer_loop = iteration_outer_loop +1 ;
    end



    iteration_inner_loop = iteration_inner_loop + 1;

end

% look at the best overall
best_alpha = max(final_round_alpha);
index = find(final_round_alpha == best_alpha);
best_povm = final_round_povm(index(1));
best_channels = final_round_channels(index(1));

fprintf("Best of %d: %f\n", MAX_ITER_OUTER_LOOP, best_alpha);

results_per_ineq{ineq_nr} = {best_alpha, index, best_povm, best_channels};

ScenarioFilename = 'scenario_joe';
filename = strcat(ScenarioFilename,'best_arxiv_112_2626','.mat');
save(filename);
