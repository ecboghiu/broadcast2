diary arxiv1112_2626

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

rng('shuffle');

ins = [2,2,2];
outs = [2,2,2];

load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads table3arXiv11122626

nr_ineqs = size(bellcoeffs_cell,2);
results = {};
best_alpha = 0;
best_inner_iter = 0;
for ineq_nr=1:nr_ineqs
    bellcoeffs = bellcoeffs_cell{ineq_nr};
    localboundNS2 = local_upper_bounds(ineq_nr);
    quantumbound = table3arXiv11122626(ineq_nr, 5);
    
    localboundBroadcast = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
    if abs(localboundBroadcast-localboundNS2)>1e-6
       fprintf("\n WARNING: you shouldn't get a broadcast local bound greater than NS2 %f %f \n", localboundBroadcast, localboundNS2); 
    end
    
    fprintf("\n\n\n This is INEQUALITY nr ineq_nr=%d, localboundNS2=%f\n\n\n", ineq_nr, localboundNS2);

    final_round_alpha = []; % store the final visibility after a round of 
    final_round_povm = {};
    final_round_channels = {};

    % LOOP PARAMETERS
    MAX_ITER_OUTER_LOOP = 1;
    ALHPA_TOL_DIST_TO_POINT = 1e-3;
    DELTA_ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
    MAX_ITER_INNER_LOOP = 200;
    ALPHA_STDEV_TOL = 1;
    ALPHA_INI_ITER = 200;
    %%%


    % INITIAL VALUES LOOP
    latest_alpha_meta = 0;
    meta_iteration = 1;

    while meta_iteration <= MAX_ITER_OUTER_LOOP
        fprintf("\nRound %d of the see-saw.\n", meta_iteration);

        % store outputs from the loop
        inner_list_of_alphas = [];
        inner_list_of_povms = {};
        inner_list_of_channels = {};

        % initial values inner loop
        iteration_inner_loop = 1;
        alpha_stdev = ALPHA_STDEV_TOL-1;
        alpha = 0;

        while iteration_inner_loop <= MAX_ITER_INNER_LOOP 
alpha=0;

alpha_iteration=0;
LPstatus = 0;
%% FIND A GOOD INITIAL POINT
while (abs(alpha-0)<ALHPA_TOL_DIST_TO_POINT || abs(alpha-1)<ALHPA_TOL_DIST_TO_POINT) && (alpha_iteration < ALPHA_INI_ITER)
            POVMs = givePprojRANDgeneral(ins);
            channel = {giveChannelRAND(2,4)};
            [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyPartiallyEntangled(1-0.7, 0.1), POVMs, channel);

            p_entangled = ProbMultidimArray(final_state(NoisyPartiallyEntangled(1, 0.1), channel), POVMs);
            p_uniform   = ProbMultidimArray(final_state(NoisyPartiallyEntangled(0, 0.1), channel), POVMs);
            alpha = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
            fprintf("ini cond alpha=%f",alpha);
            if abs(alpha-0)> ALHPA_TOL_DIST_TO_POINT && abs(alpha-1)> ALHPA_TOL_DIST_TO_POINT %finalObj-localboundNS2>1e-6
                p_entangled = ProbMultidimArray(final_state(NoisyPartiallyEntangled(1, 0.1), channel), POVMs);
                p_uniform   = ProbMultidimArray(final_state(NoisyPartiallyEntangled(0, 0.1), channel), POVMs);
                alpha = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
                fprintf("\nFound good initial condition. Bell value:%f NS2Bound: %f alpha=%f\n", finalObj, localboundNS2, alpha);
                break;
            end
            fprintf("ini cond iter = %d ineq=%d/%d meta_iter=%d/%d alpha_iter=%d/%d\n", ...
                alpha_iteration, ineq_nr, nr_ineqs, meta_iteration, MAX_ITER_OUTER_LOOP, alpha_iteration, ALPHA_INI_ITER );
            alpha_iteration = alpha_iteration + 1;
end
            %% DO SEE SAW STARTING FROM THIS POINT
            ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
            MAX_ITER = 50;
            AbsDeltaAlpha = 1e6;
            iteration = 1;
            alpha_stdev = 1;
            while (AbsDeltaAlpha>ALPHA_CONVERGENCE_THRESHOLD && ...
                    iteration<MAX_ITER && ...
                    abs(alpha-0)>ALHPA_TOL_DIST_TO_POINT)
                [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(alpha), POVMs, channel);
                oldalpha = alpha;
                p_entangled = ProbMultidimArray(final_state(NoisyPartiallyEntangled(1, 0.1), channel), POVMs);
                p_uniform   = ProbMultidimArray(final_state(NoisyPartiallyEntangled(0, 0.1), channel), POVMs);
                alpha = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
                AbsDeltaAlpha = abs(alpha-oldalpha);
                fprintf('alpha=%f ineq=%d/%d meta_iter=%d/%d inner_iter=%d/%d\n', ...
                    alpha, ineq_nr, nr_ineqs, meta_iteration, MAX_ITER_OUTER_LOOP, iteration_inner_loop, MAX_ITER_INNER_LOOP);
            end
            fprintf("alpha converged to: %f in iteration: %d ineq_nr: %d\n", alpha, best_inner_iter, ineq_nr);
            if alpha-(1-0.683)>1e-8 && abs(alpha-1)>1e-3 && abs(alpha)>1e-3
               fprintf("\n\n\nFound what we're looking for. Bell value:%f NS2 Bound: %f alpha=%f\n\n\n\n", finalObj, localboundNS2,alpha);
               save('final_workspace.mat');
               error("Finishing program.");
            end
            
            if alpha > best_alpha && alpha < 0.9
               best_alpha = alpha; 
               best_inner_iter = iteration_inner_loop;
            end
            fprintf("Best alpha so far = %f\n", best_alpha);
                
            inner_list_of_alphas = [inner_list_of_alphas, alpha];
            inner_list_of_povms{iteration_inner_loop} = POVMs;
            inner_list_of_channels{iteration_inner_loop} = channel;
            iteration_inner_loop = iteration_inner_loop + 1;

        end
        final_round_alpha = [final_round_alpha, inner_list_of_alphas(end)];
        final_round_povm{length(final_round_alpha)} = inner_list_of_povms{end};
        final_round_channels{length(final_round_alpha)} = inner_list_of_channels{end};
        latest_alpha_meta = final_round_alpha(end);
        meta_iteration = meta_iteration + 1;
    end

    % look at the best overall
    best_alpha = max(final_round_alpha);
    index = find(final_round_alpha == best_alpha);
    best_povm = final_round_povm(index(1));
    best_channels = final_round_channels(index(1));

    fprintf("Best of %d: %f\n", MAX_ITER_OUTER_LOOP, best_alpha);

    results{ineq_nr} = {best_alpha, index, best_povm, best_channels};
    
 
end
ScenarioFilename = 'scenario';
filename = strcat(ScenarioFilename,'best_arxiv_112_2626','.mat');
save(filename);
