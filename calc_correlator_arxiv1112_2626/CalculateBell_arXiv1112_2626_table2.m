mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

% Fix the seed for debugging.
%rng(1234); % this seems to saturate the local bound
rng('default');

ins = [2,2,2];
outs = [2,2,2];

% % % % A1 = sym('A1');
% % % % A2 = sym('A2');
% % % % B1 = sym('B1');
% % % % B2 = sym('B2');
% % % % C1 = sym('C1');
% % % % C2 = sym('C2');
% % % % % From the table description 'table2_arXiv1112.2626.txt'
% % % % % The order of the coefficients is the following :
% % % % %     constant, <A0>, <A1>, <B0>, <A0 B0>, <A1 B0>, <B1>, <A0 B1>, <A1 B1>,
% % % % %     <C0>, <A0 C0>, <A1 C0>, <B0 C0>, <A0 B0 C0>, <A1 B0 C0>, <B1 C0>, <A0 B1 C0>, <A1 B1 C0>,
% % % % %     <C1>, <A0 C1>, <A1 C1>, <B0 C1>, <A0 B0 C1>, <A1 B0 C1>, <B1 C1>, <A0 B1 C1>, <A1 B1 C1>
% % % % % so that the first inequality translates to
% % % % %     -1 - <A1> - <B1> - <A1 B1> - <C1> - <A1 C1> - <B1 C1> - <A1 B1 C1> <= 0
% % % % % Therefore we create a vector
% % % % 
% % % % rand_input = 1;% randi(ins(1));
% % % % Pid_A = sym(strcat('Pi_1_',string(rand_input),'_1'))+sym(strcat('Pi_1_',string(rand_input),'_2'));
% % % % rand_input = 1;% randi(ins(2));
% % % % Pid_B = sym(strcat('Pi_2_',string(rand_input),'_1'))+sym(strcat('Pi_2_',string(rand_input),'_2'));
% % % % rand_input = 1;% randi(ins(3));
% % % % Pid_C = sym(strcat('Pi_3_',string(rand_input),'_1'))+sym(strcat('Pi_3_',string(rand_input),'_2'));
% % % % 
% % % % correlator_vec_tuple = [A1*Pid_B*Pid_C, A2*Pid_B*Pid_C, Pid_A*B1*Pid_C, A1*B1*Pid_C, A2*B1*Pid_C, Pid_A*B2*Pid_C, A1*B2*Pid_C, A2*B2*Pid_C, ...
% % % %                   Pid_A*Pid_B*C1, A1*Pid_B*C1, A2*Pid_B*C1, Pid_A*B1*C1, A1*B1*C1, A1*B1*C1, Pid_A*B2*C1, A1*B2*C1, A2*B2*C1, ...
% % % %                   Pid_A*Pid_B*C2, A1*Pid_B*C2, A2*Pid_B*C2, Pid_A*B1*C2, A1*B1*C2, A2*B1*C2, Pid_A*B2*C2, A1*B2*C2, A2*B2*C2];             
% % % %               
load('table2_arXiv1112_2626.mat');
nr_inequalities = size(table2arXiv11122626,1);
inequalities = [];
local_upper_bounds = [];
% % % % 
% % % % Pi_1_1_1 = sym('Pi_1_1_1');
% % % % Pi_1_1_2 = sym('Pi_1_1_2');
% % % % Pi_1_2_1 = sym('Pi_1_2_1');
% % % % Pi_1_2_2 = sym('Pi_1_2_2');
% % % % Pi_2_1_1 = sym('Pi_2_1_1');
% % % % Pi_2_1_2 = sym('Pi_2_1_2');
% % % % Pi_2_2_1 = sym('Pi_2_2_1');
% % % % Pi_2_2_2 = sym('Pi_2_2_2');
% % % % Pi_3_1_1 = sym('Pi_3_1_1');
% % % % Pi_3_1_2 = sym('Pi_3_1_2');
% % % % Pi_3_2_1 = sym('Pi_3_2_1');
% % % % Pi_3_2_2 = sym('Pi_3_2_2');
% % % % all_vars = [A1, A2, B1, B2, C1, C2];
% % % % substitution_rules = [Pi_1_1_1-Pi_1_1_2, ...
% % % %                      Pi_1_2_1-Pi_1_2_2, ...
% % % %                      Pi_2_1_1-Pi_2_1_2, ...
% % % %                      Pi_2_2_1-Pi_2_2_2, ...
% % % %                      Pi_3_1_1-Pi_3_1_2, ...
% % % %                      Pi_3_2_1-Pi_3_2_2];

corr_to_prob_mat = fromCorrToBellMat(ins,outs);
prob_to_corr_mat = fromBellToCorrMat(ins,outs);

%fprintf('WARNING:REALLY SLOW\n');
bellcoeffs_cell = {};
for i=1:nr_inequalities
    %disp(i);
    local_upper_bounds = [local_upper_bounds, -table2arXiv11122626(i,1)];
% %     ineq = dot( table2arXiv11122626(i, 2:end), correlator_vec_tuple);
% %     ineq = subs(ineq, all_vars, substitution_rules);   
% %     ineq = expand(ineq);
% %     ineq_prob = FromProjectorNotationToProbability(ineq,ins,outs);
% %     bellcoeffs_ineq = GetBellCoeffsFromProbSymIneq(ineq_prob,ins,outs);
    bellcoeffs_cell{i} = reshape(corr_to_prob_mat*table2arXiv11122626(i, 2:end)',[ins,outs]);
    
    %% Check that my function "fromCorrToBellMat" behaves well with the inverse I wrote "fromBellToCorrMat"
    
    b_prob = bellcoeffs_cell{i};
    b_corr = table2arXiv11122626(i,2:end); % 2:end instead of 1:end because we ignore the local bound
    
    b_corr_wmat = prob_to_corr_mat * b_prob(:);
    b_prob_wmat = corr_to_prob_mat * b_corr(:);
    
    assert(norm(b_corr_wmat(2:end)-b_corr(:)) < 1e-12, "Functions are not good inverses of one another");
    assert(norm(b_prob_wmat(:)-b_prob(:)) < 1e-12, "Functions are not good inverses of one another");
end




% Don't use this it doesn't work correctly.
%LBoundI3322 = ClassicalOptInequalityNoBroadcast(bellcoeffsI3322)

save('bellcoeffs_arxiv1112_2626.mat','bellcoeffs_cell','local_upper_bounds','ins','outs');