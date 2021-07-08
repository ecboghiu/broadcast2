clear all

diary arxiv1112_2626.txt

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
addpath(newdir);
addpath(newdir2);


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
results_per_ineq = cell(NR_OF_INEQS, 3);
results_within_loop = cell(1,3);

tic
for ineq_nr=2:NR_OF_INEQS
    %% Set seed randomly so that we can reproduce the statistics within one big loop
    rng('shuffle','twister');
    fprintf("rng info:\n\n");
    disp(rng);
    
    bellcoeffs = bellcoeffs_cell{ineq_nr};
    localboundNS2 = local_upper_bounds(ineq_nr);
    quantumbound = table3arXiv11122626(ineq_nr, 5);
    
    fprintf("\n\n\tInequality number = %d\n\n\n", ineq_nr);
    
    %% Loop parameters
    ALPHA_INI_ITER = 100; % How many initial conditions to try for the initial point
    MAX_ITER_BIG_LOOP = 10; % 4How many times to try the whole process
    MAX_ITER_VIS_OPT_LOOP = 50; % given good initial condition, loop for trying to optimize the visibility
    
    INITIAL_VISIBILITY = 0;  % !! I'm using the convention of (1-p)psi+p*id instead of p*psi+(1-p)*id
    CONST_CHI = 0.15;  % For the partially entangled states
    IDENTITY_PLACEMENT = 'A'; % If using rho_A \otimes Id/2 (IDENTITY_PLACEMENT = 'B') or Id/2 \otimes rho_B (IDENTITY_PLACEMENT = 'B') 
    STATE_SETTINGS = struct('name','werner');
    
    
    % Solve ineq (23) from https://arxiv.org/pdf/1510.06721.pdf
    COS2 = (cos(2*CONST_CHI))^2;
    roots23 = roots([COS2,-2*COS2,0,2,-1]);
    roots23 = roots23(abs(imag(roots23))<1e-8); % discard imaginary roots;
    roots23 = roots23(abs(roots23)>=0);
    roots23 = roots23(abs(roots23)<=1); % discard p outside [0,1]
    roots23 = max(roots23); % just in case there is more than one, but there shouldnt be
    p = roots23 - 0.1; % check that ineq (23) is satisfied in the intervat (0,roots23)
    assert(dot([COS2,-2*COS2,0,2,-1],[p^4,p^3,p^2,p^1,1]) <= 0, "Something bad happened");
    % for all p >
    UNSTEERABILITY_THRESHOLD = 1 - roots23; % we're using opposite convention for noise
    
    ALHPA_TOL_DIST_TO_POINT = 1e-4; % Util for excluding "bad" points such 
                                    % as visibility=0 or visibility=1 (if 
                                    % its zero, then most likely its an infeasibility)
    ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
    VISIBILITY_CONVERGENCE_THRESHOLD = 1e-6;
    
    % Suppose we have a certain visibility eta in the quantum state
    % With this noisy quantum state we get a visibility nu for a bell
    % inequality
    % With the definition (1-eta)*rho+eta*id for the visibility
    % we can absorb the visibility of the bell inequality into the quantum
    % state as follows:
    % eta' = eta + nu - eta*nu
    % We can update eta -> eta' and then try to optimize the bell
    % inequality again over channels and measurements. Doing this might
    % make us fall into a local optimum. For this, when absorbing nu into
    % eta, we try to fall short a bit to help explore the landscape more.
    % This parameter controls the "effective inequality visibility", that
    % is what percentage of nu we absorb into eta.
    % eta' = eta + nu*DELTA_STATE_VIS_PROP - eta*nu*DELTA_STATE_VIS_PROP
    DELTA_STATE_VIS_PROP = 0.75;  
    
    % Calculate the best per inequality
    best_visibility = 0;
    best_visibility_iter = 1;
    best_POVMS = givePprojRANDgeneral(ins);
    best_channel = {giveChannelRAND(2,4)};
    best_visisibility_ineq_nr = 1;
    results_within_loop = cell(1);
    
    iteration_big_loop = 1;
    while iteration_big_loop <= MAX_ITER_BIG_LOOP 
        
        
        alpha=0;
        alpha_iteration=1;
        LPstatus = 0;
        
        state_visibility = INITIAL_VISIBILITY;
        %% FIND A GOOD INITIAL POINT
        fprintf("\nNext iteration of big loop.\nLooking for a good initial condition.\n");
        while (abs(alpha-0)<ALHPA_TOL_DIST_TO_POINT || abs(alpha-1)<ALHPA_TOL_DIST_TO_POINT) && (alpha_iteration <= ALPHA_INI_ITER)
            %fprintf("Initial condition search iter: %d\n", alpha_iteration);
            POVMs = givePprojRANDgeneral(ins);  % Random projective measurements
            channel = {giveChannelRAND(2,4)};  % Random isometry
            
            %[POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyState(visibility, STATE_SETTINGS), POVMs, channel);
            [POVMs,finalObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, NoisyState(state_visibility, STATE_SETTINGS), POVMs, channel, ins, outs);
            
            p_entangled = ProbMultidimArray(final_state(NoisyState(INITIAL_VISIBILITY, STATE_SETTINGS), channel), POVMs, ins, outs);
            p_uniform   = ProbMultidimArray(final_state(NoisyState(1, STATE_SETTINGS), channel), POVMs, ins, outs);
            [alpha, LPstatus]= visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
            
            if LPstatus == 0
                aux1 = bellcoeffs.*p_entangled; aux2 = bellcoeffs.*p_uniform;
                fprintf("Ineq.nr.:%d/%d Big.loop:%d/%d Ini.iter.:%d/%d Vis of ineq given random ini = %f. vis_state=%f Bell·p_ent=%f Bell·p_unif=%f bound_l=%f bound_q=%f \n", ...
                    ineq_nr, NR_OF_INEQS, iteration_big_loop, MAX_ITER_BIG_LOOP, alpha_iteration, ALPHA_INI_ITER, alpha, state_visibility, sum(aux1(:)), sum(aux2(:)), localboundNS2, quantumbound);
                if abs(alpha-0)> ALHPA_TOL_DIST_TO_POINT && abs(alpha-1)> ALHPA_TOL_DIST_TO_POINT
                    fprintf("   Found good initial condition. Bell value:%f NS2Bound:%f ineq_visibility= %f state_vis=%f\n", finalObj, localboundNS2, alpha, state_visibility);
                    break;
                end
            else
               warning("Visibility LP not solver correctly (LPstatus=%d). Trying another initial condition.", LPstatus);
            end
            alpha_iteration = alpha_iteration + 1;
        end
        if alpha_iteration > ALPHA_INI_ITER
           warning("Couldn't find a good initial point before hitting the maximum number of iterations. Going to a different inequality.")
           break;
        else
            fprintf("   Given this point, optimize the visibility of the inequality.\n");
            %% OPTIMIZE OVER THE VISIBILITY
            iter_vis_opt_loop = 1;
            old_visibility = -1;
            AbsDeltaStateVis = 1e6;
            while iter_vis_opt_loop <= MAX_ITER_VIS_OPT_LOOP && AbsDeltaStateVis > VISIBILITY_CONVERGENCE_THRESHOLD
                [POVMs,finalObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, NoisyState(state_visibility, STATE_SETTINGS), POVMs, channel, ins, outs);
                channel = cleanChannel(channel, 2, 4);
                %checkPOVMsAreGood(POVMs,ins,outs);
                %checkThatChannelIsGood(channel, 2, 4);
                
                p_entangled = ProbMultidimArray(final_state(NoisyState(state_visibility, STATE_SETTINGS), channel), POVMs, ins, outs);
                p_uniform   = ProbMultidimArray(final_state(NoisyState(1,                STATE_SETTINGS), channel), POVMs, ins, outs);
                [ineq_visibility, LPstatus]= visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
             
                % If LP infeasible, or visibilities are ridiculous, add
                % a little bit of noise to the POVMS and channel and try
                % again. If this still happens reduce slightly the noise in
                % the quantum state. If it still happens discard round.
                if abs(ineq_visibility-1)<ALHPA_TOL_DIST_TO_POINT || ineq_visibility > 1 || ineq_visibility < 0 || LPstatus ~= 0
                    % try to add a bit of white noise to see if the
                    % situation improves
                    noise = 1e-6;
                    warning("Getting ridiculous visibilities (vis=%f LPstatus=%d).\nTrying to add white noise = %g to the channel and measurements and reducing noise in the quantum state by %g.\n", ineq_visibility, LPstatus, noise, noise);
                    
                    state_visibility_lessvisibility = state_visibility - noise;
                    if state_visibility_lessvisibility < 0
                       state_visibility_lessvisibility = state_visibility_lessvisibility + noise; % if we go below 0, just don't do the substraction 
                    end
                    
                    channel_noise = (1-noise) * channel + noise * eye(size(channel,1));
                    POVMs_noise = POVMs;
                    for party=1:nrparties
                        for x=1:ins(party)
                            for a=1:outs(party)
                                POVMs_noise{party}{x}{a} = (1-noise) * POVMs_noise{party}{x}{a} + noise * eye(size(POVMs_noise{party}{x}{a},1));
                            end
                        end
                    end
                    p_entangled = ProbMultidimArray(final_state(NoisyState(state_visibility, STATE_SETTINGS), channel_noise), POVMs_noise, ins, outs);
                    p_uniform   = ProbMultidimArray(final_state(NoisyState(1, STATE_SETTINGS), channel_noise), POVMs_noise, ins, outs);
                    [ineq_visibility, LPstatus] = visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
                    
                    warning("After adding white noise and reduing state_visibility we get vis=%f LPstatus=%d", ineq_visibility, LPstatus);
                    
                    if abs(ineq_visibility-1)<ALHPA_TOL_DIST_TO_POINT || ineq_visibility > 1 || ineq_visibility < 0 || LPstatus ~= 0
                        warning("Adding noise didn't work, restarting with a different initial condition.");
                        break;
                    end
                end

% % %                 % Then we go back a bit and decrease the noise and try again
% % %                 if state_visibility + ineq_visibility*DELTA_STATE_VIS_PROP - state_visibility*ineq_visibility*DELTA_STATE_VIS_PROP > 0 && abs(ineq_visibility-0) > ALHPA_TOL_DIST_TO_POINT
% % %                     % If the visibility of the ineq is different from 0
% % %                     % then we update the state visibility a bit below that
% % %                     % of the inequality visibility
                old_visibility = state_visibility;
                %fprintf("Ineq vis converged to: %f. Update state vis from %f to %f. DeltaVis=%f\n", ineq_visibility, state_visibility, state_visibility + ineq_visibility - DELTA_STATE_VIS,ineq_visibility - DELTA_STATE_VIS);
                state_visibility = state_visibility + ineq_visibility*DELTA_STATE_VIS_PROP - state_visibility*ineq_visibility*DELTA_STATE_VIS_PROP; % remember I use the opposite convention for visibility
                AbsDeltaStateVis = abs(state_visibility-old_visibility);
% % %                 else
% % %                     % if we are here we are at the beginning of the loop
% % %                     % where doing the DELTA_VIS CHANGE might actually give
% % %                     % us negative visibilities
% % %                     % if so we just updat to ineq_visibility and not
% % %                     % slightly below, since going slightly below gives us
% % %                     % negative values
% % %                     old_visibility = state_visibility;
% % %                     %fprintf("Ineq vis converged to: %f. Update state vis from %f to %f. DeltaVis=%f\n", ineq_visibility, state_visibility, state_visibility + ineq_visibility,ineq_visibility);
% % %                     state_visibility = state_visibility + ineq_visibility - state_visibility*ineq_visibility; % remember I use the opposite convention for visibility
% % %                     AbsDeltaStateVis = abs(state_visibility-old_visibility);
% % %                 end

                
                fprintf("Ineq.nr.:%d/%d Big.loop.:%d/%d Vis_iter=%d Ineq vis after optimizing: %f with vis_state=%f. Updating vis_state from %f to %f. (Bellvalue=%f bound_l=%f bound_q=%f).\n", ...
                        ineq_nr, NR_OF_INEQS, iteration_big_loop, MAX_ITER_BIG_LOOP, iter_vis_opt_loop, ineq_visibility, old_visibility, old_visibility, state_visibility, finalObj, localboundNS2, quantumbound);

                iter_vis_opt_loop = iter_vis_opt_loop + 1;
            end
            if iter_vis_opt_loop <= MAX_ITER_VIS_OPT_LOOP
                %[~,finalObj,~] = optimizer_seesaw(optimizer_objects, bellcoeffs, NoisyState(0, STATE_SETTINGS), POVMs, channel, ins, outs);
                fprintf("State visibility converged to %.10g +- %g.  Starting a new iteration.\n", state_visibility, VISIBILITY_CONVERGENCE_THRESHOLD);
            else
                warning("State visibility didn't converge within the maximum iterations. Trying another initial condition.\n"); 
            end
        end
        if state_visibility > best_visibility
            best_visibility = state_visibility;
            best_visibility_iter = iteration_big_loop;
            best_POVMS = POVMs;
            best_channel = channel;
            best_visisibility_ineq_nr = ineq_nr;
            
            if strcmp(STATE_SETTINGS.name, 'partially_entangled')
                if STATE_SETTINGS.IDENTITY_PLACEMENT == 'B'
                    fprintf("\n\tNew best visibility found! p=%f (UNS A->B, Steerable B->A for p>=%f)\n", state_visibility, UNSTEERABILITY_THRESHOLD);
                elseif STATE_SETTINGS.IDENTITY_PLACEMENT == 'A'
                    fprintf("\n\tNew best visibility found! p=%f (UNS B->A, Steerable A->B for p>=%f)\n", state_visibility, UNSTEERABILITY_THRESHOLD);
                end
            elseif strcmp(STATE_SETTINGS.name, 'werner')
                fprintf("\n\tNew best visibility found! p=%f\n", state_visibility);
            end
                
            fprintf("Channel:\n");
            disp(best_channel);
            fprintf("Channel spectrum=\n");
            disp(eig(best_channel).');
            fprintf("Measurements:\n");
            for party=1:nrparties
                for x=1:ins(party)
                    obs_x = POVMs{party}{x}{1}-POVMs{party}{x}{2};
                    bloch = BlochComponents(obs_x);
                    bloch = num2cell(bloch(2:4));
                    [azimuth,elevation,r] = cart2sph(bloch{:});
                    azimuth = azimuth*180/pi;
                    elevation = elevation*180/pi;
                    fprintf("Party: %d, input:%d, obs: (azimuth[º],elevation[º],r)=(%f,%f,%f)\n", party, x, azimuth, elevation, r);
                end
            end
           
        end
        results_within_loop{iteration_big_loop, 1} = state_visibility;
        
        iteration_big_loop = iteration_big_loop + 1;
    end
    all_vis_for_ineq = [results_within_loop{:,1}];
    promedio = mean(all_vis_for_ineq);
    desviacion = sqrt(var(all_vis_for_ineq));
    
    results_per_ineq{ineq_nr, 1} = best_visibility;
    results_per_ineq{ineq_nr, 2} = best_channel;
    results_per_ineq{ineq_nr, 3} = best_POVMS;
    results_per_ineq{ineq_nr, 4} = promedio;
    results_per_ineq{ineq_nr, 5} = desviacion; 
    results_per_ineq{ineq_nr, 6} = [results_within_loop{:,1}]; 
    
    all_vis = [results_per_ineq{:,1}];
    wheremax = all_vis == max(all_vis);
    auxints = 1:length(all_vis);
    position_max = auxints(wheremax); % this might be degenerate i just take the first
    value_max = all_vis(wheremax);
    
    
    
    fprintf("\n\t Starting new inequality. Current best is state_visibility=%f for inequality %d. Average: %g +- %g (1 sigma).\n\n", value_max(1), position_max(1), promedio, desviacion);
    

    % Save to file
    fprintf("Saving current workspace to file.\n");
    ScenarioFilename = 'scenario';
    filename = strcat(ScenarioFilename,'best_arxiv_1112_2626','.mat');
    save(filename);
    
end
toc


plot([results_per_ineq{:,1}])

function [POVMs,newObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, state_AB, POVMs, channel, ins, outs)

optimizer_ch = optimizer_objects{1};
optimizer_a = optimizer_objects{2};
optimizer_b = optimizer_objects{3};
optimizer_c = optimizer_objects{4};

%% Loop parameters
MAX_NR_ITERATIONS = 500;
CONVERGENCE_TOL = 1e-6;
%%

deltaObj = 1e6;
oldObj = -1e6;
newObj = 0;

iteration = 1;
while deltaObj > CONVERGENCE_TOL && iteration <= MAX_NR_ITERATIONS 
    oldObj = newObj;
    %% Channel
    belloperator = give_Bell_operator(bellcoeffs, POVMs, ins, outs);
    state = state_AB;
    ia_state = give_ia_state(state);

    %output = optimizer_ch([{ia_state}, {belloperator}]);
    output = optimizer_ch({ia_state,belloperator});
    channel = output{1};
    
%     newObj = output{2};

    %% Do the POVMs
    output_state = final_state(state, channel);
    ia_output_state = give_ia_state(output_state);

    %% Alice
    partial_products_for_a = give_partial_products(POVMs, bellcoeffs, 1, ins, outs);
    output = optimizer_a([{ia_output_state}, partial_products_for_a(:)']);
%     oldObj = newObj;
%     newObj = output{1};
%     deltaObj = abs(newObj-oldObj);
    povm_alice = reshape({output{2:end}},[ins(1),outs(1)]);
    for x=1:ins(1)
        for a=1:outs(1)
            POVMs{1}{x}{a} = povm_alice{x,a};  % Update the POVMs
        end
    end

    %% Bob
    partial_products_for_b = give_partial_products(POVMs, bellcoeffs, 2, ins, outs);
    output = optimizer_b([{ia_output_state}, partial_products_for_b(:)']);
%     oldObj = newObj;
%     newObj = output{1};
%     deltaObj = abs(newObj-oldObj);
    povm_bob = reshape({output{2:end}},[ins(2),outs(2)]);
    for y=1:ins(2)
        for b=1:outs(2)
            POVMs{2}{y}{b} = povm_bob{y,b};  % Update the POVMs
        end
    end

    %% Charlie
    partial_products_for_c = give_partial_products(POVMs, bellcoeffs, 3, ins, outs);
    output = optimizer_c([{ia_output_state}, partial_products_for_c(:)']);
    
    povm_charlie = reshape({output{2:end}},[ins(3),outs(3)]);
    for z=1:ins(3)
        for c=1:outs(3)
            POVMs{3}{z}{c} = povm_charlie{z,c};  % Update the POVMs
        end
    end
    
    newObj = output{1};
    deltaObj = abs(newObj-oldObj);

    iteration = iteration + 1;
end
if iteration > MAX_NR_ITERATIONS
   warning("Hit maximum number of iterations (=%d) without conv. with tol=%g. Returning last value as converged.", MAX_NR_ITERATIONS, CONVERGENCE_TOL); 
end

end
