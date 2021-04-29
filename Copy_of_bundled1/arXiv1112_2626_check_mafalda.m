clear all

diary arxiv1112_2626.txt

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
addpath(newdir);
addpath(newdir2);

%% Set seed
rng('shuffle','twister');
fprintf("rng info:\n\n");
disp(rng);

%% Scenario settings
load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads 'table3arXiv11122626'

ins = [2,2,2];
outs = [2,2,2];
nrparties = length(ins);

%% Calculate the 'optimizer' objects
optimizer_ch = optimizer_channel(ins, outs, 2);
optimizer_a = optimizer_povm_party(1, ins, outs);
optimizer_b = optimizer_povm_party(2, ins, outs);
optimizer_c = optimizer_povm_party(3, ins, outs);
optimizer_objects = {optimizer_ch, optimizer_a, optimizer_b, optimizer_c};

%% 
NR_OF_INEQS = size(bellcoeffs_cell,2);
results_per_ineq = cell(1,NR_OF_INEQS);
for ineq_nr=2:10
    bellcoeffs = bellcoeffs_cell{ineq_nr};
    localboundNS2 = local_upper_bounds(ineq_nr);
    quantumbound = table3arXiv11122626(ineq_nr, 5);
    
    fprintf("\n\n\tInequality number = %d\n\n\n", ineq_nr);
    
    %% Loop parameters
    ALPHA_INI_ITER = 10; % How many initial conditions to try for the initial point
    MAX_ITER_BIG_LOOP = 10; % How many times to try the whole process
    MAX_ITER_VIS_OPT_LOOP = 50; % given good initial condition, loop for trying to optimize the visibility
    
    INITIAL_VISIBILITY = 0;  % !! I'm using the convention of (1-p)psi+p*id instead of p*psi+(1-p)*id
    CONST_CHI = 0.01;  % For the partially entangled states
    
    ALHPA_TOL_DIST_TO_POINT = 1e-4; % Util for excluding "bad" points such 
                                    % as visibility=0 or visibility=1 (if 
                                    % its zero, then most likely its an infeasibility)
    ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
    VISIBILITY_CONVERGENCE_THRESHOLD = 1e-3;
    
    DELTA_STATE_VIS = 0.0;  % when updating the state visibility to that 
                             % of the ineq visibility, we go slightly 
                             % below and this controls that

    iteration_big_loop = 1;
    while iteration_big_loop <= MAX_ITER_BIG_LOOP 
        alpha=0;
        alpha_iteration=1;
        LPstatus = 0;

        state_visibility = INITIAL_VISIBILITY;
        %% FIND A GOOD INITIAL POINT
        fprintf("Looking for a good initial condition.\n");
        while (abs(alpha-0)<ALHPA_TOL_DIST_TO_POINT || abs(alpha-1)<ALHPA_TOL_DIST_TO_POINT) && (alpha_iteration <= ALPHA_INI_ITER)
            fprintf("Initial condition search iter: %d\n", alpha_iteration);
            POVMs = givePprojRANDgeneral(ins);  % Random projective measurements
            channel = {giveChannelRAND(2,4)};  % Random isometry
            
            %[POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyPartiallyEntangled(visibility, CONST_CHI), POVMs, channel);
            [POVMs,finalObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, NoisyPartiallyEntangled(state_visibility, CONST_CHI), POVMs, channel, ins, outs);
            
            p_entangled = ProbMultidimArray(final_state(NoisyPartiallyEntangled(INITIAL_VISIBILITY, CONST_CHI), channel), POVMs);
            p_uniform   = ProbMultidimArray(final_state(NoisyPartiallyEntangled(1, CONST_CHI), channel), POVMs);
            alpha = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
            
            aux1 = bellcoeffs.*p_entangled; aux2 = bellcoeffs.*p_uniform;
            fprintf("Vis of ineq given random ini = %f. vis_state=%f Bell·p_ent=%f Bell·p_unif=%f bound_l=%f bound_q=%f Ineq.nr.:%d/%d Ini.iter.:%d/%d\n", ...
                alpha, state_visibility, sum(aux1(:)), sum(aux2(:)), localboundNS2, quantumbound, ineq_nr, NR_OF_INEQS, alpha_iteration, ALPHA_INI_ITER);
            if abs(alpha-0)> ALHPA_TOL_DIST_TO_POINT && abs(alpha-1)> ALHPA_TOL_DIST_TO_POINT
                fprintf("\nFound good initial condition. Bell value:%f NS2Bound:%f alpha= %f state_vis=%f\n", finalObj, localboundNS2, alpha, state_visibility);
                break;
            end
            alpha_iteration = alpha_iteration + 1;
        end
        if alpha_iteration > ALPHA_INI_ITER
           warning("Couldn't find a good initial point before hitting the maximum number of iterations. Going to a different inequality.")
           break;
        else
            fprintf("Given this point, optimize the visibility of the inequality.\n\n");
            %% OPTIMIZE OVER THE VISIBILITY
            iter_vis_opt_loop = 1;
            old_visibility = -1;
            AbsDeltaStateVis = 1e6;
            while iter_vis_opt_loop <= MAX_ITER_VIS_OPT_LOOP && AbsDeltaStateVis > VISIBILITY_CONVERGENCE_THRESHOLD
                [POVMs,finalObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, NoisyPartiallyEntangled(state_visibility, CONST_CHI), POVMs, channel, ins, outs);
                p_entangled = ProbMultidimArray(final_state(NoisyPartiallyEntangled(state_visibility, CONST_CHI), channel), POVMs);
                p_uniform   = ProbMultidimArray(final_state(NoisyPartiallyEntangled(1, CONST_CHI), channel), POVMs);
                ineq_visibility = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
                if abs(ineq_visibility-1)<ALHPA_TOL_DIST_TO_POINT
                    % try to add a bit of white noise to see if the
                    % situation improves
                    noise = 1e-6;
                    channel_noise = (1-noise) * channel + noise * channel;
                    for party=1:nrparties
                        for x=1:ins(party)
                            for a=1:outs(party)
                                POVMs{party}{x}{a} = (1-noise) * POVMs{party}{x}{a} + noise * POVMs{party}{x}{a};
                            end
                        end
                    end
                    p_entangled = ProbMultidimArray(final_state(NoisyPartiallyEntangled(state_visibility, CONST_CHI), channel), POVMs);
                    p_uniform   = ProbMultidimArray(final_state(NoisyPartiallyEntangled(1, CONST_CHI), channel), POVMs);
                    ineq_visibility = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
                    if abs(ineq_visibility-1)<ALHPA_TOL_DIST_TO_POINT
                        warning("Getting visibility=1 which means that probably both probabilities are inside the local set. Restarting with a different initial condition..");
                        break;
                    end
                end

                % Then we go back a bit and decrease the noise and try again
                if state_visibility + ineq_visibility - DELTA_STATE_VIS > 0 && abs(ineq_visibility-0) > ALHPA_TOL_DIST_TO_POINT
                    % If the visibility of the ineq is different from 0
                    % then we update the state visibility a bit below that
                    % of the inequality visibility
                    old_visibility = state_visibility;
                    %fprintf("Ineq vis converged to: %f. Update state vis from %f to %f. DeltaVis=%f\n", ineq_visibility, state_visibility, state_visibility + ineq_visibility - DELTA_STATE_VIS,ineq_visibility - DELTA_STATE_VIS);
                    state_visibility = state_visibility + ineq_visibility - state_visibility*ineq_visibility - DELTA_STATE_VIS; % remember I use the opposite convention for visibility
                    AbsDeltaStateVis = abs(state_visibility-old_visibility);
                else
                    % if we are here we are at the beginning of the loop
                    % where doing the DELTA_VIS CHANGE might actually give
                    % us negative visibilities
                    % if so we just updat to ineq_visibility and not
                    % slightly below, since going slightly below gives us
                    % negative values
                    old_visibility = state_visibility;
                    %fprintf("Ineq vis converged to: %f. Update state vis from %f to %f. DeltaVis=%f\n", ineq_visibility, state_visibility, state_visibility + ineq_visibility,ineq_visibility);
                    state_visibility = state_visibility + ineq_visibility - state_visibility*ineq_visibility; % remember I use the opposite convention for visibility
                    AbsDeltaStateVis = abs(state_visibility-old_visibility);
                end
                if state_visibility > 0.3
                   disp(1);
                end
                
                fprintf("Ineq vis after optimizing: %f with vis_state=%f. Updating vis_state from %f to %f. (Bellvalue=%f bound_l=%f bound_q=%f). Vis_iter=%d. Ineq %d/%d\n", ...
                 ineq_visibility, old_visibility, old_visibility, state_visibility, finalObj, localboundNS2, quantumbound, iter_vis_opt_loop, ineq_nr, NR_OF_INEQS);

                iter_vis_opt_loop = iter_vis_opt_loop + 1;
            end
            aux1 = bellcoeffs .* ProbMultidimArray(final_state(NoisyPartiallyEntangled(0, CONST_CHI), channel), POVMs);
            fprintf("State visibility converged to %f. Ineq. value with no state noise: %f (max quantum: %f). Starting a new iteration.\n\n", state_visibility, sum(aux1(:)), quantumbound);
        end
        iteration_big_loop = iteration_big_loop + 1;
    end

end

% Save to file
ScenarioFilename = 'scenario';
filename = strcat(ScenarioFilename,'best_arxiv_112_2626','.mat');
save(filename);

function [POVMs,newObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, state_AB, POVMs, channel, ins, outs)

optimizer_ch = optimizer_objects{1};
optimizer_a = optimizer_objects{2};
optimizer_b = optimizer_objects{3};
optimizer_c = optimizer_objects{4};

%% Loop parameters
MAX_NR_ITERATIONS = 200;
CONVERGENCE_TOL = 1e-7;
%%

deltaObj = 1e6;
oldObj = -1e6;
newObj = 0;

iteration = 1;
while deltaObj > CONVERGENCE_TOL && iteration <= MAX_NR_ITERATIONS   
    %% Channel
    belloperator = give_Bell_operator(bellcoeffs, POVMs, ins, outs);
    state = state_AB;
    ia_state = give_ia_state(state);

    output = optimizer_ch([{ia_state}, {belloperator}]);
    channel = output{1};
    oldObj = newObj;
    newObj = output{2};

    %% Do the POVMs
    output_state = final_state(state, channel);
    ia_output_state = give_ia_state(output_state);

    %% Alice
    partial_products_for_a = give_partial_products(POVMs, bellcoeffs, 1, ins, outs);
    output = optimizer_a([{ia_output_state}, partial_products_for_a(:)']);
    oldObj = newObj;
    newObj = output{1};
    deltaObj = abs(newObj-oldObj);
    povm_alice = reshape({output{2:end}},[ins(1),outs(1)]);
    for x=1:ins(1)
        for a=1:outs(1)
            POVMs{1}{x}{a} = povm_alice{x,a};  % Update the POVMs
        end
    end

    %% Bob
    partial_products_for_b = give_partial_products(POVMs, bellcoeffs, 2, ins, outs);
    output = optimizer_b([{ia_output_state}, partial_products_for_b(:)']);
    oldObj = newObj;
    newObj = output{1};
    deltaObj = abs(newObj-oldObj);
    povm_bob = reshape({output{2:end}},[ins(2),outs(2)]);
    for y=1:ins(2)
        for b=1:outs(2)
            POVMs{2}{y}{b} = povm_bob{y,b};  % Update the POVMs
        end
    end

    %% Charlie
    partial_products_for_c = give_partial_products(POVMs, bellcoeffs, 3, ins, outs);
    output = optimizer_c([{ia_output_state}, partial_products_for_c(:)']);
    oldObj = newObj;
    newObj = output{1};
    deltaObj = abs(newObj-oldObj);
    povm_charlie = reshape({output{2:end}},[ins(3),outs(3)]);
    for z=1:ins(3)
        for c=1:outs(3)
            POVMs{3}{z}{c} = povm_charlie{z,c};  % Update the POVMs
        end
    end

iteration = iteration + 1;
if iteration == MAX_NR_ITERATIONS
   warning("Hit maximum number of iterations. Returning last value as converged."); 
end
end

end
