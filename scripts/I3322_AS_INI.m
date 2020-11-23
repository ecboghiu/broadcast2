diary mymatlabdiary

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

MAX_ITER_META = 3;

%%%%%%%%%
% Fix the seed for debugging.
%rng(1234); % this seems to saturate the local bound
rng('default');


ScenarioFilename = 'I3322broadcast';
load(strcat('bellcoeffs',ScenarioFilename,'.mat')); % this loads 'bellcoeffs' 'ins' 'outs'
LBoundOurEBI = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
fprintf('Local bound: %f\n' ,LBoundOurEBI);

final_round_alpha = [];
final_round_povm = {};
final_round_channels = {};

latest_alpha_meta = 0;

meta_iteration = 1;
while meta_iteration < MAX_ITER_META && latest_alpha_meta < 1-0.577
    fprintf("\nRound %d of the see-saw.\n", meta_iteration);
    % While the visibility is close to 1 or above 1, keep changing the 
    % initial conditions until we get a behaviour outside the broadcast-local
    % polytope.
    alpha = 0;
    ALHPA_TOL_DIST_TO_1 = 1e-3;
    LPstatus = 0;
    [localbound, LPstatus] = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
    fprintf("Local bound of the bell inequality: %f\n", localbound);
    
    while abs(alpha-0)<ALHPA_TOL_DIST_TO_1 || LPstatus~=0
        % first 'kick start' the maximization of the inequality
        spv0=-1;
        spv1=0;
        % the point of this lopp is to get AT LEAST get a starting point
        % that isn't worse than just the maximally mixed state
        while abs(spv0-spv1)<1e-3
            %POVMs = givePprojRAND2(ins,outs);
            POVMs = givePprojRANDgeneral(ins);
            %POVMs = givePprojRANDmaxEBI();
            %POVMs = givePprojDET();
            %channel = RandomSuperoperator([2,4]);
            channel = {giveChannelRAND(2,4)};
            %channel = {give_Joe_U()};
            %channel = giveChannelAddsIdentity(2,2,"right");
            
            spv0 = evaluate_bell_ineq(bellcoeffs, 0, final_state(NoisyWernerState(0), channel), POVMs);
            spv1 = evaluate_bell_ineq(bellcoeffs, 0, final_state(NoisyWernerState(1), channel), POVMs);
            fprintf("s·p(v=0)=%f s·p(v=1)=%f localbound=%f\n", spv0, spv1, localbound);
        end
        p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs);
        p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs);
        alpha = visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform);
        
        % once we have a point s.t. b·p_ent > b·p_uni we do quantum
        % optimization
        [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(alpha), POVMs, channel);
        assert(checkPOVMsAreGood(POVMs,ins,outs),'Problem with POVMs');
        assert(checkThatChannelIsGood(channel, 2, 4), 'Problem with the channel');
        
        % calculate the visibility of this, if the visibility is not
        % greater than 0 then start again
        p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs);
        p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs);
        alpha = visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform);
        fprintf("Local maximum: %f\n", alpha);
    end

    % Parameters for stopping conditions for the next loop.
    ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
    MAX_ITER = 1000;
    AbsDeltaAlpha = 1e6;
    iteration = 1;
    alpha_stdev = 1;

    list_of_alphas = [];
    list_of_povms = {};
    list_of_channels = {};
    while (AbsDeltaAlpha>ALPHA_CONVERGENCE_THRESHOLD && ...
            iteration<MAX_ITER && ...
            abs(alpha-0)>ALHPA_TOL_DIST_TO_1 && ...
            alpha_stdev > 1e-3)
        %[localbound, LPstatus] = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
        %fprintf("Local bound of bell inequality: %f\n", localbound);
            fprintf("s·p(v=0), s·p(v=alpha), s·p(v=1), s·p1, s·p2, s·(p1-p2): %f, %f, %f, %f, %f, %f\n", ...
        evaluate_bell_ineq(bellcoeffs, 0, final_state(NoisyWernerState(0), channel), POVMs), ...
        evaluate_bell_ineq(bellcoeffs, 0, final_state(NoisyWernerState(alpha), channel), POVMs), ...
        evaluate_bell_ineq(bellcoeffs, 0, final_state(NoisyWernerState(1), channel), POVMs),...
        sum(bellcoeffs .* (p_entangled),'all'), ...
        sum(bellcoeffs .* (p_uniform),'all'), ...
        sum(bellcoeffs .* (p_entangled-p_uniform),'all'));
        localbound = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
        fprintf("Visibility check: %f\n", visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform));
        [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(alpha), POVMs, channel);
%         POVMs = newPOVMS;
%         channel = newChannel;
        assert(checkPOVMsAreGood(POVMs,ins,outs),'Problem with POVMs');
        assert(checkThatChannelIsGood(channel, 2, 4), 'Problem with the channel');
        

        %outputcritvis = criticalvisibility_std(Pproj, channel, ins, outs);
        oldalpha = alpha;
        p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs);
        p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs);
        if norm(p_entangled(:) - p_uniform(:)) < 1e-6
            fprintf("Distributions too equal to each other. Starting again.\n");
           break; % the LP would not work so it's best to start over
        end
        [newalpha, newbellcoeffs, LPstatus, dual_alpha] = BroadcastLP(p_entangled, p_uniform);
        if LPstatus ~= 0
            %disp('LP not solved correctly');
            fprintf("LP not solved correctly. Trying another set of initial points.\n");
            break;
            %error('Check why infeasible problem.'); 
            list_of_alphas = [list_of_alphas, 0];
            list_of_povms{iteration} = 0;
            list_of_channels{iteration} = 0;
        else
            bellcoeffs = newbellcoeffs;
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
            list_of_channels{iteration} = channel;
            iteration = iteration + 1;
                
            AbsDeltaAlpha = abs(alpha-oldalpha);

            fprintf("Visibility after optimizing over POVMs and channels: %f\n", alpha);
        end
        %alpha = outputcritvis{1} + 0.01;
    end
    if ~isempty(list_of_alphas)
        final_round_alpha = [final_round_alpha, list_of_alphas(end)];
        final_round_povm{length(final_round_alpha)} = list_of_povms{end};
        final_round_channels{length(final_round_alpha)} = list_of_channels{end};
        latest_alpha_meta = final_round_alpha(end);
    end
    meta_iteration = meta_iteration + 1;
end

best_alpha = max(final_round_alpha);
index = find(final_round_alpha == best_alpha);
best_povm = final_round_povm(index);
best_channels = final_round_channels(index);

fprintf("Best of %d: %f 1-%f=%f\n", MAX_ITER_OUTER_LOOP, best_alpha,best_alpha,1-best_alpha);

save('matlabworkspace_asini.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%save 'data_optimal0.577Standard'
% 
% LBoundOurEBI = ClassicalOptInequality2(bellcoeffs,ins,outs);
% fprintf('LBoundOurEBI: %f\n', LBoundOurEBI);
% 
% inputstate = final_state(ini_state(alpha), channel);
% mixingstate = eye(size(inputstate))/size(inputstate,1);
% val=evaluate_bell_ineq(bellcoeffs, 0, inputstate, Pproj, ins, outs);
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
% ineqconstant+evaluate_bell_ineq(bellcoeffsagain, ineqconstant, final_state(ini_state(0.577), channel), Pproj, ins, outs)
% 
% LBoundOurEBI2 = ClassicalOptInequality2(bellcoeffsagain,ins,outs)
% fprintf("local bounds before and after: %f %f\n", LBoundOurEBI, LBoundOurEBI2);
