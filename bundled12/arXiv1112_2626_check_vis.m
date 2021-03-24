diary arxiv1112_2626.txt

clear all

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
addpath(newdir);
addpath(newdir2);

rng('shuffle');

ins = [2,2,2];
outs = [2,2,2];

load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads 'table3arXiv11122626'

best_alpha = 0;
best_inner_iter = 0;
best_visibility = 0;

NR_OF_INEQS = size(bellcoeffs_cell,2);
results_per_ineq = cell(1,NR_OF_INEQS);
for ineq_nr=120:130
    bellcoeffs = bellcoeffs_cell{ineq_nr};
    localboundNS2 = local_upper_bounds(ineq_nr);
    quantumbound = table3arXiv11122626(ineq_nr, 5);
    
    fprintf("\n\n\n Inequality number = %d\n\n\n", ineq_nr);

    %% Sanity check with my own code
    localboundBroadcast = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
    if abs(localboundBroadcast-localboundNS2)>1e-6
       warning("You shouldn't get a broadcast local bound greater than NS2.  BroadcastL = %f NS2L = %f \n", localboundBroadcast, localboundNS2); 
    end
    
    %% Loop parameters
    ALPHA_INI_ITER = 50; % How many initial conditions to try for the initial point
    MAX_ITER_INNER_LOOP = 200; % How many times to try the whole process
    ALHPA_TOL_DIST_TO_POINT = 1e-3; % Util for excluding "bad" points such as visibility=0 or visibility=1 (if its zero, then most likely its an infeasibility)
    ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
    MAX_ITER_INNER = 50;
    MAX_ITER_OUTER = 50;
    INITIAL_VISIBILITY = 0;  % !! I'm using the convention of (1-p)psi+p*id instead of p*psi+(1-p)*id
    CONST_CHI = 0.1;
    DELTA_STATE_VIS = 0.01;
    
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

    iteration_inner_loop = 1;
    while iteration_inner_loop <= MAX_ITER_INNER_LOOP 
        alpha=0;
        alpha_iteration=0;
        LPstatus = 0;

        visibility = INITIAL_VISIBILITY;
        %% FIND A GOOD INITIAL POINT
        while (abs(alpha-0)<ALHPA_TOL_DIST_TO_POINT || abs(alpha-1)<ALHPA_TOL_DIST_TO_POINT) && (alpha_iteration < ALPHA_INI_ITER)
            POVMs = givePprojRANDgeneral(ins);  % Random projective measurements
            channel = {giveChannelRAND(2,4)};  % Random isometry
            [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyPartiallyEntangled(visibility, CONST_CHI), POVMs, channel);
            p_entangled = ProbMultidimArray(final_state(NoisyPartiallyEntangled(INITIAL_VISIBILITY, CONST_CHI), channel), POVMs);
            p_uniform   = ProbMultidimArray(final_state(NoisyPartiallyEntangled(1, CONST_CHI), channel), POVMs);
            alpha = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
            
            aux1 = bellcoeffs.*p_entangled; aux2 = bellcoeffs.*p_uniform;
            fprintf("Vis of ineq given random ini = %f. vis_state=%f Bell·p_ent=%f Bell·p_unif=%f bound_l=%f bound_q=%f Ineq.nr.:%d/%d Ini.iter.:%d/%d\n", ...
                alpha, visibility, sum(aux1(:)), sum(aux2(:)), localboundNS2, quantumbound, ineq_nr, NR_OF_INEQS, alpha_iteration, ALPHA_INI_ITER);
            if abs(alpha-0)> ALHPA_TOL_DIST_TO_POINT && abs(alpha-1)> ALHPA_TOL_DIST_TO_POINT && alpha < 0.1
                fprintf("\nFound good initial condition. Bell value:%f NS2Bound: %f alpha=%f state_vis=%f\n", finalObj, localboundNS2, alpha, visibility);
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
                [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyPartiallyEntangled(visibility, CONST_CHI), POVMs, channel);
                p_entangled = ProbMultidimArray(final_state(NoisyPartiallyEntangled(visibility, CONST_CHI), channel), POVMs);
                p_uniform   = ProbMultidimArray(final_state(NoisyPartiallyEntangled(1, CONST_CHI), channel), POVMs);
                alpha = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
                AbsDeltaAlpha = abs(alpha-oldalpha);
                %fprintf("Visibility of bell inequality after optimizing = %f. Bellvalue=%f localbound=%f quantumboud=%f ineq. nr.:%d/%d\n", alpha, finalObj, localboundNS2, quantumbound, ineq_nr, NR_OF_INEQS);
                aux1 = bellcoeffs.*p_entangled; aux2 = bellcoeffs.*p_uniform;
                fprintf("Visibility of ineq %d/%d after optimizing = %f. vis_state=%f  Bellvalue=%f Bell·p_ent=%f Bell·p_unif=%f bound_l=%f bound_q=%f\n", ...
                ineq_nr, NR_OF_INEQS, alpha, visibility, finalObj, sum(aux1(:)), sum(aux2(:)), localboundNS2, quantumbound);
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
               fprintf("Ineq vis converged to 0, reached local optimum. Restart with different initial conditions.\n\n");
               break;
            end
            % Then we go back a bit and decrease the noise and try again
            old_visibility = visibility;
            fprintf("Ineq vis converged to: %f. Update state vis from %f to %f. DeltaVis=%f\n", alpha, visibility, visibility + alpha - DELTA_STATE_VIS, + alpha - DELTA_STATE_VIS);
            visibility = visibility + alpha - DELTA_STATE_VIS; % remember I use the opposite convention for visibility
            AbsDeltaStateVis = abs(visibility-old_visibility);
            iteration_outer_loop = iteration_outer_loop +1 ;
            
            if visibility > best_visibility
                best_visibility = visibility;
            end
            fprintf("Current best state visibility p = %f\n", best_visibility);
        end

        final_round_alpha = inner_list_of_alphas(end);
        final_round_povm{iteration_inner_loop} = inner_list_of_povms{end};
        final_round_channels{iteration_inner_loop} = inner_list_of_channels{end};
        
        iteration_inner_loop = iteration_inner_loop + 1;

    end

    % look at the best overall
    best_alpha = max(final_round_alpha);
    index = find(final_round_alpha == best_alpha);
    best_povm = final_round_povm(index(1));
    best_channels = final_round_channels(index(1));

    fprintf("Best of %d: %f\n", MAX_ITER_OUTER_LOOP, best_alpha);

    results_per_ineq{ineq_nr} = {best_alpha, index, best_povm, best_channels};
end
ScenarioFilename = 'scenario';
filename = strcat(ScenarioFilename,'best_arxiv_112_2626','.mat');
save(filename);
