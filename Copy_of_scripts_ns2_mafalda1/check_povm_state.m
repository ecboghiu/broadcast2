load('chiresults.mat');
% 'chiresults' is a cell that containts:
% chiresults{:,1} the chi values
% chiresults{:,2} the p_crit for which tr Pi^xyz_abc rho_POVM(p_crit,chi) = p(abc|xyz) has for Bell ineq 16 as score the local bound 
% chiresults{:,3} the p_threshold for one way steerability. we want p_crit < p_threshold
% chiresults{:,4} p_crit again because I'm dumb
% chiresults{:,5} the channel in choi form
% chiresults{:,6} the measurements. the projectors onto outputs are called as nested cells {party}{input}{output}

% Tune which point in the graph to look at
CHI_IDX = 2;

chi         = chiresults{CHI_IDX,1};
channel     = chiresults{CHI_IDX,5};
channel = cleanChannel(channel, 2, 4);
disp(min(eig(channel)))
disp(norm(PartialTrace(channel, 2, [2,4])-eye(2)))

POVMs       = chiresults{CHI_IDX,6};
for p=1:length(ins)
   for x=1:ins(p)
      for a=1:outs(p)
          if ~IsPSD(POVMs{p}{x}{a})
              %warning("POVM not positive! min(eig)=%g  p,x,a=%d,%d,%d\n", min(eig(POVMs{p}{x}{a})), p, x, a);
              POVMs{p}{x}{a} = (1-eta)*POVMs{p}{x}{a} + eta*eye(2);
          end
      end
      summ = 0;
      for a=1:outs(p)
          summ = summ + POVMs{p}{x}{a};
      end
      POVMs{p}{x}{2} = POVMs{p}{x}{a} + eye(2) - summ;
%       for a=1:outs(p)
%           POVMs{p}{x}{a} = POVMs{p}{x}{a}*2/trace(summ) ;
%       end
%           if norm(summ-eye(2),'fro')>1e-12
%           warning("POVMs don't summ up to 1! dist_to_eye(2)=%g\n", norm(summ-eye(2)));
%           for a=1:outs(p)
%               %POVMs{p}{x}{a} = POVMs{p}{x}{a}/trace(summ);
%           end
%           end
   end
end
fprintf("new")
for p=1:length(ins)
   for x=1:ins(p)
       fprintf(" p x: %d %d \n",p,x);
       summ = 0;
      for a=1:outs(p)
           summ = summ + POVMs{p}{x}{a};
            if min(eig(POVMs{p}{x}{a})) < 0
                fprintf("flag")
                v=abs(min(eig(POVMs{p}{x}{a})));
               disp(min(eig(POVMs{p}{x}{a}))) 
               POVMs{p}{x}{a} = (1-v)*POVMs{p}{x}{a}+v*eye(2);
               fprintf("new eig: %f\n", min(eig(POVMs{p}{x}{a})));
            end
          if ~IsPSD(POVMs{p}{x}{a})
              warning("POVM not positive! min(eig)=%g  p,x,a=%d,%d,%d\n", min(eig(POVMs{p}{x}{a})), p, x, a);
          end
      end
        fprintf("dist to eye: %g\n", (norm(summ-eye(2))))
   end
end



p_threshold = chiresults{CHI_IDX,3};

% Uncomment to see that the non-unitarity of the channel is important (for
% these exact measurements)
kr = KrausOperators(channel,[2,4]);
thrash = zeros(4,2);

%kr{1} = thrash;
%kr{2} = thrash;
%kr{3} = thrash;
%kr{4} = thrash;
%kr{5} = thrash;
%kr{6} = thrash;
%kr{7} = thrash;
channel = ChoiMatrix(kr);


%% Scenario settings
load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads 'table3arXiv11122626'

ins = [2,2,2];
outs = [2,2,2];
nrparties = length(ins);

ineq_nr = 16;

bellcoeffs = bellcoeffs_cell{ineq_nr};
localboundNS2 = local_upper_bounds(ineq_nr);
quantumbound = table3arXiv11122626(ineq_nr, 5);

% WARNING: I use the opposite convention for visibiltiies. For example for 
% the isotropic state it p=0 is no white noise, and p=1 is max white noise.
p_entangled = ProbMultidimArray(final_state(PartiallyEntangledPOVM(0, chi, 'A'), channel), POVMs, ins, outs);
p_uniform   = ProbMultidimArray(final_state(PartiallyEntangledPOVM(1, chi, 'A'), channel), POVMs, ins, outs);
p_uniform2   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs, ins, outs);

b1 = bellcoeffs.*p_entangled;
b2 = bellcoeffs.*p_uniform;
b3 = bellcoeffs.*p_uniform2;
b1 = sum(b1(:));
b2 = sum(b2(:));
b3 = sum(b3(:));

% NOTE: In the following function I do the following I find the p for
% which the Bell inequality has the NS2 bound value.
% find p s.t. (1-p) b·p_entangled + p b·p_uniform == localboundNS2
% and I call this p_critical
% then I absorpb the p into a p for the POVM state:
% (1-p) b·p_entangled + p b·p_uniform = b·p(p_effective)
% I checked that the distribution p(vis,chi) that results from the POVM_state(vis,chi) is linear in
% vis and as such the p_effective == p_critical
[p_crit, ~]= visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);

fprintf("p_crit=%f and p_onewaysteerability=%f (we want p_crit > p_onewaysteerability) b1 b2 b3 = %f %f %f\n", p_crit, p_threshold, b1, b2, b3);

%% Uncomment to plot the observables
% Check for different CHI_IDX!
% BlochAllObs(POVMs)

%% Plot p_crit vs chi
% figure
% hold on
% plot([chiresults{:,1}],[chiresults{:,2}])
% plot([chiresults{:,1}],[chiresults{:,3}])
% hold off

function state = PartiallyEntangledPOVM(p,xi,identity_placement)
    e1 = [1;0];
    e2 = [0;1];

    e1e1 = Tensor(e1,e1);
    e1e2 = Tensor(e1,e2);
    e2e1 = Tensor(e2,e1);
    e2e2 = Tensor(e2,e2);
    
%     e1e1 = [1;0;0;0];
%     e2e2 = [0;0;0;1];
    
    psi_xi = cos(xi) * e1e1 + sin(xi) * e2e2;
    rho_xi = psi_xi * psi_xi'/ (psi_xi' * psi_xi);
    
    ket00 = zeros(2,2);
    ket00(1,1) = 1;
    
    if identity_placement == 'A'
        rhoB = PartialTrace(rho_xi, 1, [2,2]);
        state = (1-p) * rho_xi + p * kron(eye(2)/2, rhoB);
        state_A =  PartialTrace(state, 2, [2,2]);
        state = 0.5 * state + 0.5 * kron(state_A, ket00);
    elseif identity_placement == 'B'
        rhoA = PartialTrace(rho_xi, 2, [2,2]);
        state = (1-p) * rho_xi + p * kron(rhoA, eye(2)/2);
        state_B =  PartialTrace(state, 1, [2,2]);
        state = 0.5 * state + 0.5 * kron(ket00, state_B);
    else
       error("Invalid 'identity_placement' input."); 
    end
    
    assert(IsPSD(state),"Not positive!");
    assert(trace(state)-1<1e-12,"Not normalized!");
end

function state = final_state(inistate,channel)
    % TODO put dimA, dimB, dimB1, dimB2 as function inputs only works for
    % 2,2,2 scenario
    dimA=2;
    dimB=2;
    dimB1=2;
    dimB2=2;
    
    biggerstate = kron( inistate.', eye(dimA*dimB1*dimB2) );
    
    Phi = auxPHI(dimA);
    biggerchannel = kron(Phi*Phi',ChoiMatrix(channel));
    swapop=Tensor(eye(dimA),SwapOperator(2),eye(dimB1*dimB2));
    biggerchannel = swapop * biggerchannel * swapop';
    
    % TODO Rewrite this to use only 1 dimA, inefficient for more outputs
    state = PartialTrace(biggerchannel*biggerstate, [1,2], [dimA,dimB,dimA,dimB1,dimB2]);
end

function PHI = auxPHI(dim)
    PHI = zeros(dim*dim,1);
    basis=eye(dim);
    for i=1:dim
       PHI = PHI + kron(basis(:,i),basis(:,i));
    end
end

function probability_ndarray = ProbMultidimArray(state,povms, inputs_per_party, outputs_per_party)
    nrparties = length(inputs_per_party);
    dims = num2cell([inputs_per_party,outputs_per_party]);
    probability_ndarray = zeros(dims{:});
    aux = [inputs_per_party, outputs_per_party];
    allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
    for slice=1:size(allinputoutputcombinations,1)
        ins = num2cell(allinputoutputcombinations(slice,1:nrparties));
        outs = num2cell(allinputoutputcombinations(slice,nrparties+1:end));
        tensor = 1;
        for p=1:nrparties
             tensor = kron(tensor,povms{p}{ins{p}}{outs{p}});
        end
        probability_ndarray(ins{:},outs{:}) = real(trace(tensor*state));
    end
end

function [visibility, LPstatus] = visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform)

%% LP way
% assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
% dims = size(bellcoeffs);
% nrparties = length(dims)/2;
% ins = dims(1:nrparties);
% outs = dims(nrparties+1:end);
% 
% alpha = sdpvar(1);
% vis_constraints = [alpha>=0, alpha<=1];
% 
% p_noisy = (1-alpha)*p_entangled + alpha*p_uniform;
% % objective = sum(bellcoeffs .* p_noisy, size(bellcoeffs));
% % aux = bellcoeffs .* p_noisy;
% % auxsize=size(aux);
% % while prod(auxsize(:)) ~= 1
% %    aux = sum(aux); 
% % end
% 
% objective = 0;
% coords = ind2subv(size(bellcoeffs), 1:prod(dims(:)));
% for idx = 1:size(coords,1)
%     coords_choice = num2cell(coords(idx,:));
%     objective = objective + bellcoeffs(coords_choice{:}) * p_noisy(coords_choice{:});
% end
% %sdisplay(objective);
% constraints = [vis_constraints, objective <= localbound];
% 
% optsol = optimize(constraints, -objective, sdpsettings('solver','mosek','verbose', 0));
% LPstatus = optsol.problem;
% if optsol.problem ~= 0
%     warning('Infeasibility in inequality visibility LP.');
%    %error('check why this LP is not working'); 
% end
% 
% visibility = value(alpha);

%% Manual way

aux1 = bellcoeffs .* p_entangled;
c1 = sum(aux1(:));

aux2 = bellcoeffs .* p_uniform;
c2 = sum(aux2(:));

comparison_tolerance = 1e-6;
if abs(c1-c2) < comparison_tolerance
    visibility = 0;
    LPstatus = 1;
    diff_p1p2 = norm(p_entangled(:) - p_uniform(:));
    if diff_p1p2 < comparison_tolerance
        warning("The probability distributions are too close. tol: %f", comparison_tolerance);
    else
        warning("The Bell values of the probability distributions are too close. tol: %f", comparison_tolerance);
    end
    return;
end

if c1 <= localbound && c2 <= localbound
    visibility = 0;
    % technically it is an infeasible problem, but for this case it's best 
    % not to raise this flag as this is a common situation. 
    % it usually gives visibility=1 or somthing like that
    LPstatus = 0; 
    
    %warning("Neither p1 nor p2 violate the Bell inequality. Returning 0.");
    return;
elseif c1 >= localbound && c2 >= localbound
    visibility = 0;
    LPstatus = 1;
    warning("Both p1 and p2 violate the Bell inequality. Returning 0.");
    return;
elseif c1 <= localbound && c2 >= localbound
    visibility = 0;
    LPstatus = 1;
    warning("The code convention is that p1 is outside the local set and p2 inside. You probably didn't optimize over channels and measurements. Returning 0.");
    return;
else
   % only thing left is c1 >= localbound and c2 <= localbound
   visibility = (c1-localbound)/(c1-c2);
   LPstatus = 0;
end

end

function subscript_vector = ind2subv(siz, linear_index)
%IND2SUBV return vector position of a linear index into an ND-array
%
% subscript_vector = ind2subv(siz, linear_index)
%
% Inputs:
%                   siz Dx1 or 1xD   size of array to index
%          linear_index Nx1 or 1xN   required elements using Fortran order
%
% Outputs:
%     subscript_vector      NxD      required elements as subscript row-vectors
%
% Matlab's ind2sub turns linear scalar indexes into the subscripts corresponding
% to the row, col, 3rd-dim, ... of the element. The subscript for each dim
% is returned as a separate argument. This version returns a single argument:
% each row is a vector containing all the subscripts for the corresponding
% linear_index.
%
% See also: subv2ind, ind2sub, sub2ind

% Iain Murray, November 2007, March 2011
% This technique appeared earlier in the following forum post:
%
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5476&objectType=File
% Date: 2004-10-04
% From: Mukhtar Ullah (mukhtar.ullah@informatik.uni-rostock.de)

[out{1:length(siz)}] = ind2sub(siz, linear_index(:));
subscript_vector = cell2mat(out);
end

function BlochAllObs(povms)
    nr_parties = size(povms,2);
    nr_inputs_per_party = [];
    nr_outputs_per_party = [];
    for party=1:nr_parties
        inps = size(povms{party},2);
        outs = size(povms{party}{1},1);
        nr_inputs_per_party = [nr_inputs_per_party, inps];
    end
    
    bloch_components = {{}};
    string_names = {{}};
    for party=1:nr_parties
       party_char = char(party+'A'-1);
       for x=1:nr_inputs_per_party(party)
           party_in = string(x);
           string_names{party}{x} = strcat(party_char, party_in);
           allbloch = BlochComponents(povms{party}{x}{1}-povms{party}{x}{2});
           bloch_components{party}{x} = allbloch(2:end); % remove id
       end
    end
    
    
    figure(1)
    hold on
    
    [Xs, Yx, Zx] = sphere(25);
    mySphere = surf( Xs, Yx, Zx );
    axis equal
    shading interp
    xlabel('x')
    ylabel('y')
    zlabel('z')
    mySphere.FaceAlpha = 0.25;
    
    line([-1 1], [0 0], [0 0])
    line([0 0], [-1 1], [0 0])
    line([0 0], [0 0], [-1 1])
    
    text(0, 0, 1.1, '$\left| 0 \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')
    text(1.1, 0, 0, '$\left| + \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')
    text(-1.1, 0, 0, '$\left| - \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')    
    text(0, 0, -1.1, '$\left| 1 \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')
    

    view([60 12])
    party_colors = {'r','b','g'};
    past_vectors = [[0;0;0]];
    for party=1:nr_parties
        for x=1:nr_inputs_per_party(party)
            p1 = [0 0 0];
            p2 = bloch_components{party}{x} ;
            dp = p2-p1;
%             past_vectors = [past_vectors, p2'];
%             flag = false;
%             while flag == false
%                 p2 = p2 + [0.1,0.1,0.1];
%                 for vec=1:size(past_vectors,2)
%                     if past_vectors(:,vec) == p2
%                     flag = false;
%                     break;
%                     end
%                 end
%                 disp(past_vectors);
%             end
            text(p2(1)+0.1, p2(2)+0.1, p2(3)+0.1, string_names{party}{x}, 'FontSize', 14, 'HorizontalAlignment', 'Center')
            
            quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),'LineWidth',3,'DisplayName',string_names{party}{x},'Color',party_colors{party})

        end
    end
    
    legend()
    
end

function out = BlochComponents(obs)
    sig0 = eye(2);
    sig1 = [0+0i,1+0i;1+0i,0+0i];
    sig2 = [0+0i,-1i;1i,0+0i];
    sig3 = [1+0i,0+0i;0+0i,-1+0i] ;
    
    r0 = trace(obs' * sig0);
    r1 = trace(obs' * sig1);
    r2 = trace(obs' * sig2);
    r3 = trace(obs' * sig3);
    
    out = [r0, r1, r2, r3];
end


