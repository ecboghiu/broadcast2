load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads table3arXiv11122626

ineq_nr = 99;
bellcoeffs_ref = bellcoeffs_cell{ineq_nr};
localbound = local_upper_bounds(ineq_nr);

assert(mod(length(size(bellcoeffs_ref)),2)==0,"There should be as many inputs as outputs.");
dims = size(bellcoeffs_ref);
nrparties = length(dims)/2;
ins = dims(1:nrparties);
outs = dims(nrparties+1:end);

POVMs = givePprojRANDgeneral(ins);
% refactoredPOVMs = {};
% for party = 1:nrparties
%    for input = 1:ins(party)
%       for output = 1:outs(party)
%          refactoredPOVMs{party,input,output,1} = real(POVMs{party}{input}{output}); 
%          refactoredPOVMs{party,input,output,2} = imag(POVMs{party}{input}{output}); 
%       end
%    end
% end
channel = {giveChannelRAND(2,4)};

%