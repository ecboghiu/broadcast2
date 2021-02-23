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
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) + p(11|00)
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
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) + p(11|00)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(1,1,1,1)=1;
prob(1,1,1,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*A1 == prob(:)), "Incorrect mat");

%% B1
prob = zeros(dimcell{:});
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) + p(11|00)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(1,1,1,1)=1;
prob(1,1,1,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*A1 == prob(:)), "Incorrect mat");

%% A1*B1
prob = zeros(dimcell{:});
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) + p(11|00)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(1,1,1,1)=1;
prob(1,1,1,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*A1 == prob(:)), "Incorrect mat");

%% A2*B1
prob = zeros(dimcell{:});
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) + p(11|00)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(1,1,1,1)=1;
prob(1,1,1,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*A1 == prob(:)), "Incorrect mat");

%% B2
prob = zeros(dimcell{:});
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) + p(11|00)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(1,1,1,1)=1;
prob(1,1,1,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*A1 == prob(:)), "Incorrect mat");

%% A1*B2
prob = zeros(dimcell{:});
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) + p(11|00)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(1,1,1,1)=1;
prob(1,1,1,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*A1 == prob(:)), "Incorrect mat");

%% B2*B2
prob = zeros(dimcell{:});
% A1 = \sum_ab (-1)^a p(ab|00) = p(00|00) + p(01|00) - p(10|00) + p(11|00)
% We choose y=1 for convention. Since it's NS it doesn't matter.
% A1 is input x=1 for x in {1,2}. So actually in the above notation
% we have inputs in {0,1} and outputs in {0,1} as well. So take this into
% account.
prob(1,1,1,1)=1;
prob(1,1,1,2)=1;
prob(1,1,2,1)=-1;
prob(1,1,2,2)=-1;
assert(all(mat_chsh*A1 == prob(:)), "Incorrect mat");