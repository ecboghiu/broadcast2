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

A1 = sym('A1');
A2 = sym('A2');
B1 = sym('B1');
B2 = sym('B2');
C1 = sym('C1');
C2 = sym('C2');
% From the table description 'table2_arXiv1112.2626.txt'
% The order of the coefficients is the following :
%     constant, <A0>, <A1>, <B0>, <A0 B0>, <A1 B0>, <B1>, <A0 B1>, <A1 B1>,
%     <C0>, <A0 C0>, <A1 C0>, <B0 C0>, <A0 B0 C0>, <A1 B0 C0>, <B1 C0>, <A0 B1 C0>, <A1 B1 C0>,
%     <C1>, <A0 C1>, <A1 C1>, <B0 C1>, <A0 B0 C1>, <A1 B0 C1>, <B1 C1>, <A0 B1 C1>, <A1 B1 C1>
% so that the first inequality translates to
%     -1 - <A1> - <B1> - <A1 B1> - <C1> - <A1 C1> - <B1 C1> - <A1 B1 C1> <= 0
% Therefore we create a vector

rand_input = 1;% randi(ins(1));
Pid_A = sym(strcat('Pi_1_',string(rand_input),'_1'))+sym(strcat('Pi_1_',string(rand_input),'_2'));
rand_input = 1;% randi(ins(2));
Pid_B = sym(strcat('Pi_2_',string(rand_input),'_1'))+sym(strcat('Pi_2_',string(rand_input),'_2'));
rand_input = 1;% randi(ins(3));
Pid_C = sym(strcat('Pi_3_',string(rand_input),'_1'))+sym(strcat('Pi_3_',string(rand_input),'_2'));

correlator_vec_tuple = [A1, A2, B1, A1*B1, A2*B1, B2, A1*B2, A2*B2, ...
                  C1, A1*C1, A2*C1, B1*C1, A1*B1*C1, A1*B1*C1, B2*C1, A1*B2*C1, A2*B2*C1, ...
                  C2, A1*C2, A2*C2, B1*C2, A1*B1*C2, A2*B1*C2, B2*C2, A1*B2*C2, A2*B2*C2];             
              
load('table2_arXiv1112_2626.mat');
nr_inequalities = size(table2arXiv11122626,1);
inequalities = [];
local_upper_bounds = [];

Pi_1_1_1 = sym('Pi_1_1_1');
Pi_1_1_2 = sym('Pi_1_1_2');
Pi_1_2_1 = sym('Pi_1_2_1');
Pi_1_2_2 = sym('Pi_1_2_2');
Pi_2_1_1 = sym('Pi_2_1_1');
Pi_2_1_2 = sym('Pi_2_1_2');
Pi_2_2_1 = sym('Pi_2_2_1');
Pi_2_2_2 = sym('Pi_2_2_2');
Pi_3_1_1 = sym('Pi_3_1_1');
Pi_3_1_2 = sym('Pi_3_1_2');
Pi_3_2_1 = sym('Pi_3_2_1');
Pi_3_2_2 = sym('Pi_3_2_2');
all_vars = [A1, A2, B1, B2, C1, C2];
substitution_rules = [Pi_1_1_1-Pi_1_1_2, ...
                     Pi_1_2_1-Pi_1_2_2, ...
                     Pi_2_1_1-Pi_2_1_2, ...
                     Pi_2_2_1-Pi_2_2_2, ...
                     Pi_3_1_1-Pi_3_1_2, ...
                     Pi_3_2_1-Pi_3_2_2];


fprintf('WARNING:REALLY SLOW TODO change to using the definitions of the correlators\n');
bellcoeffs_cell = {};
for i=1:nr_inequalities
    disp(i);
    local_upper_bounds = [local_upper_bounds, -table2arXiv11122626(i,1)];
    ineq = dot( table2arXiv11122626(i, 2:end), correlator_vec_tuple);
    ineq = subs(ineq, all_vars, substitution_rules);   
    ineq = expand(ineq);
    ineq_prob = FromProjectorNotationToProbability(ineq,ins,outs);
    bellcoeffs_ineq = GetBellCoeffsFromProbSymIneq(ineq_prob,ins,outs);
    bellcoeffs_cell{i} = bellcoeffs_ineq;
end

save('bellcoeffs_arxiv1112_2626.mat','bellcoeffs_cell','local_upper_bounds','ins','outs');