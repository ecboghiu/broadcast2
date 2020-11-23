rng('default');

C1 = sym('C1');
C2 = sym('C2');
B1 = sym('B1');
B2 = sym('B2');
B3 = sym('B3');
B4 = sym('B4');
A1 = sym('A1');
A2 = sym('A2');
A3 = sym('A3');
A4 = sym('A4');

ins = [3,2,2];
outs = [2,2,2];
PAPER = A1*B1*C1 + A1*B2*C2 + A2*B2*C2 - A2*B1*C1 + A1*B1*C2 + A1*B2*C1 + A2*B1*C2 - A2*B2*C1 - 2*A3*B1 + 2*A3*B1;
PAPERinprob = ToProbabilityNotationIneqSym(PAPER,ins,outs);
PAPERbellcoeffs = GetBellCoeffsFromProbSymIneq(PAPERinprob,ins,outs);
LBoundPAPER_NUM = ClassicalOptInequality2(PAPERbellcoeffs,ins,outs)

%BigU = kron(eye(2),U);
ChoiU = ChoiMatrix({give_Joe_U()});
Pproj = givePprojDET();
channel = ChoiU;
inputstate = final_state(ini_state(0.577),channel);
alpha = VisibilityOfBellIneq(PAPERbellcoeffs, LBoundPAPER_NUM, Pproj, inputstate, ins, outs);
fprintf('LBoundPAPER_NUM alpha = %f %f\n', LBoundPAPER_NUM, alpha);
% RoundOut = CompleteRoundFixedBell(EBIbellcoeffs, LBoundOurEBI, ins, outs);

%out = MetaCompleteRoundFixedBellCorrIneqSym(OurEBI, ins, outs);
%out = MetaCompleteRoundFixedBellCorrIneqSym(OurEBI2, ins, outs);
% 
% out = MetaCompleteRoundFixedBellCorrIneqSym(PAPER, ins, outs);
% outalpha = out{1};
% outchannel = out{2};
% outPproj = out{3};
% outbell = out{4};

val=evaluate_bell_ineq(PAPERbellcoeffs, inputstate, Pproj, ins, outs);
prob = giveProbNDArray(inputstate, Pproj, ins, outs);
checkThatProbSumsToOne(prob,ins,outs);
disp([val,alpha*4*sqrt(3)]);
% disp(ClassicalOptInequality2(outbell,ins,outs));