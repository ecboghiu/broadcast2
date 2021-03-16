mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);


ins = [2,2];
outs = [2,2];

% A1*1 A2*1 1*B1 A1*B1 A2*B1 1*B2 A1*B2 A2*B2 
mat_chsh = fromCorrToBellMat(ins,outs);

A1 = [1,0,0,0,0,0,0,0]';
A2 = [0,1,0,0,0,0,0,0]';
B1 = [0,0,1,0,0,0,0,0]';
A1B1 = [0,0,0,1,0,0,0,0]';
A2B1 = [0,0,0,0,1,0,0,0]';
B2 = [0,0,0,0,0,1,0,0]';
A1B2 = [0,0,0,0,0,0,1,0]';
A2B2 = [0,0,0,0,0,0,0,1]';

dimcell = num2cell([ins,outs]);
prob = zeros(dimcell{:}); % prob(x,y,a,b)

%% A1
prob = zeros(dimcell{:});
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) - p(11|00)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(1,1,1,1)=1;
prob(1,1,1,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*A1 == prob(:)), "Incorrect mat");

%% A2
prob = zeros(dimcell{:});
% A2 = \sum_ab (-1)^a p(ab|10) = p(00|10) + p(01|10) - p(10|10) - p(11|10)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(2,1,1,1)=1;
prob(2,1,1,2)=1;
prob(2,1,2,1)=-1;
prob(2,1,2,2)=-1;
assert(all(mat_chsh*A2 == prob(:)), "Incorrect mat");

%% B1
prob = zeros(dimcell{:});
% B1 = \sum_ab (-1)^b p(ab|00) = p(00|00) + p(10|00) - p(00|00) - p(11|00)
prob(1,1,1,1)=1;
prob(1,1,2,1)=1;
prob(1,1,1,2)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*B1 == prob(:)), "Incorrect mat");

%% A1*B1
prob = zeros(dimcell{:});
% \sum_ab (-1)^a+b p(ab|00) = p(00|00) + p(11|00) - p(10|00) - p(01|00)
prob(1,1,1,1)=1;
prob(1,1,2,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,1,2)=-1;
assert(all(mat_chsh*A1B1 == prob(:)), "Incorrect mat");

%% A2*B1
prob = zeros(dimcell{:});
prob(2,1,1,1)=1;
prob(2,1,2,2)=1;
prob(2,1,2,1)=-1;
prob(2,1,1,2)=-1;
assert(all(mat_chsh*A2B1 == prob(:)), "Incorrect mat");

%% B2
prob = zeros(dimcell{:});
prob(1,2,1,1)=1;
prob(1,2,2,1)=1;
prob(1,2,1,2)=-1;
prob(1,2,2,2)=-1;
assert(all(mat_chsh*B2 == prob(:)), "Incorrect mat");

%% A1*B2
prob = zeros(dimcell{:});
prob(1,2,1,1)=1;
prob(1,2,2,2)=1;
prob(1,2,2,1)=-1;
prob(1,2,1,2)=-1;
assert(all(mat_chsh*A1B2 == prob(:)), "Incorrect mat");

%% A2*B2
prob = zeros(dimcell{:});
prob(2,2,1,1)=1;
prob(2,2,2,2)=1;
prob(2,2,2,1)=-1;
prob(2,2,1,2)=-1;
assert(all(mat_chsh*A2B2 == prob(:)), "Incorrect mat");

%% A1*B2*C2
ins = [2,2,2];
outs = [2,2,2];

dimcell = num2cell([ins,outs]);

% correlator_vec_tuple = [A1, A2, B1, A1*B1, A2*B1, B2, A1*B2, A2*B2, ...
%                   C1, A1*C1, A2*C1, B1*C1, A1*B1*C1, A1*B1*C1, B2*C1, A1*B2*C1, A2*B2*C1, ...
%                   C2, A1*C2, A2*C2, B1*C2, A1*B1*C2, A2*B1*C2, B2*C2, A1*B2*C2, A2*B2*C2];  
mat_3party = fromCorrToBellMat(ins,outs);
A1B2C1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, ...
          0, 0, 0, 0, 0, 0, 1, 0, 0, ...
          0, 0, 0, 0, 0, 0, 0, 0]';
      
% A1B2C1 = \sum_abc (-1)^a+b+c p(abc|010) = 
% = (depends how many output=1 there are, for one 1, it's -, for two 1, its
% + and for three 1 it's minus
% 000 - 001 - 010 - 100 + 011 + 101 + 110 - 111 where for eg 010 is short for
% p(010|010)
prob = zeros(dimcell{:});
prob(1,2,1,1,1,1)=1;
prob(1,2,1,1,1,2)=-1;
prob(1,2,1,1,2,1)=-1;
prob(1,2,1,2,1,1)=-1;
prob(1,2,1,1,2,2)=1;
prob(1,2,1,2,1,2)=1;
prob(1,2,1,2,2,1)=1;
prob(1,2,1,2,2,2)=-1;
assert(all(mat_3party*A1B2C1 == prob(:)), "Incorrect mat");