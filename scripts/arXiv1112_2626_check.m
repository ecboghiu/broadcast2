diary arxiv1112_2626
yalmip("clear");

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

rng('shuffle');

load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads table3arXiv11122626

nr_ineqs = size(bellcoeffs_cell,2);
results = {};
for ineq_nr=1:nr_ineqs
    bellcoeffs = bellcoeffs_cell{ineq_nr};
    localboundNS2 = local_upper_bounds(ineq_nr);
    localboundBroadcast = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
    quantumbound = table3arXiv11122626(ineq_nr, 5);
    
    if abs(localboundBroadcast-localboundNS2)>1e-6
       fprintf("\n WARNING: you shouldn't get a broadcast local bound greater than NS2 %f %f \n", localboundBroadcast, localboundNS2); 
    end
    
    fprintf("\n\n\n This is INEQUALITY nr ineq_nr=%d, localboundNS2=%f quantumbound=%f\n\n\n", ineq_nr, localboundNS2, quantumbound);

    final_round_alpha = []; % store the final visibility after a round of 
    final_round_povm = {};
    final_round_channels = {};

    % LOOP PARAMETERS
    MAX_ITER_OUTER_LOOP = 1;
    ALHPA_TOL_DIST_TO_POINT = 1e-3;
    DELTA_ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
    MAX_ITER_INNER_LOOP = 100;
    ALPHA_STDEV_TOL = 1;
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
        iteration = 1;
        alpha_stdev = ALPHA_STDEV_TOL-1;
        alpha = 0;

        while iteration < MAX_ITER_INNER_LOOP
            POVMs = givePprojRANDgeneral(ins);
            channel = {giveChannelRAND(2,4)};
            [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(1-0.8), POVMs, channel);
            if finalObj-localboundNS2>1e-6
               fprintf("\n\n\nFound what we're looking for. Bell value:%f NS2 Bound: %f\n\n\n\n", finalObj, localboundNS2);
               save('best_arxiv_112_2626.mat');
               warning("Found violation!");
               error("Finishing program.");
            end
            fprintf("Ineq = %d Iter = %d Bell value: %f NS2bound: %f and quantumbound: %f\n", ineq_nr, iteration, finalObj, localboundNS2, quantumbound);          
                
            inner_list_of_alphas = [inner_list_of_alphas, alpha];
            inner_list_of_povms{iteration} = POVMs;
            inner_list_of_channels{iteration} = channel;
            iteration = iteration + 1;

        end
        final_round_alpha = [final_round_alpha, inner_list_of_alphas(end)];
        final_round_povm{length(final_round_alpha)} = inner_list_of_povms{end};
        final_round_channels{length(final_round_alpha)} = inner_list_of_channels{end};
        latest_alpha_meta = final_round_alpha(end);
        meta_iteration = meta_iteration + 1;
        yalmip("clear");
        clear inner_list_of_alphas inner_list_of_povms inner_list_of_channels
    end

    % look at the best overall
    best_alpha = max(final_round_alpha);
    index = find(final_round_alpha == best_alpha);
    best_povm = final_round_povm(index(1));
    best_channels = final_round_channels(index(1));

    fprintf("Best of %d: %f\n", MAX_ITER_OUTER_LOOP, best_alpha);

    results{ineq_nr} = {best_alpha, index, best_povm, best_channels};
    
    clear final_round_alpha final_round_povm final_round_channels
end
ScenarioFilename = 'vis075';
filename = strcat(ScenarioFilename,'best_arxiv_112_2626','.mat');
save(filename);
