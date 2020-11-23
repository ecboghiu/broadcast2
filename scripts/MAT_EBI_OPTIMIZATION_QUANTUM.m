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

ins = [4,4,2];
outs = [2,2,2];

EBI = A1*B1 + A1*B2 - A1*B3 - A1*B4 + A2*B1 - A2*B2 + ...
      A2*B3 - A2*B4 + A3*B1 - A3*B2 - A3*B3 + A3*B4;
EBI2 = B1*A1 + B1*A2 - B1*A3 - B1*A4 + B2*A1 - B2*A2 + ...
      B2*A3 - B2*A4 + B3*A1 - B3*A2 - B3*A3 + B3*A4;
LBoundEBI = 6;
OurEBI = EBI * (C1 + C2) + LBoundEBI * A4 *(C1 - C2);
OurEBI = expand(OurEBI);
OurEBI2 = EBI2 * (C1 + C2) + LBoundEBI * B4 * (C1 - C2);
OurEBI2 = expand(OurEBI2);
OurEBIInProb = ToProbabilityNotationIneqSym(OurEBI,ins,outs);
EBIbellcoeffs = GetBellCoeffsFromProbSymIneq(OurEBIInProb,ins,outs);
%LBoundOurEBI = ClassicalOptInequality2(OurEBI,ins,outs);
LBoundOurEBI = ClassicalOptInequality2(EBIbellcoeffs,ins,outs);
fprintf('LBoundOurEBI LBoundOurEBIprime: %f %f\n' ,LBoundOurEBI,LBoundOurEBI);

% RoundOut = CompleteRoundFixedBell(EBIbellcoeffs, LBoundOurEBI, ins, outs);

%iniP_proj = givePprojRAND2();
%iniP_proj = givePprojRANDgeneral(length(ins),ins,outs);
iniP_proj = givePprojRANDmaxEBI();
%iniP_proj = givePprojDET();

%iniChannel = RandomSuperoperator([2,4]);
%iniChannel = {giveChannelRAND(2,4)};
iniChannel = {give_Joe_U()};
out = MetaCompleteRoundFixedBellCorrIneqSym(OurEBI2, iniP_proj, iniChannel, ins, outs);
%out = MetaCompleteRoundFixedBellCorrIneqSym(OurEBI2, iniP_proj, iniChannel, ins, outs);

out = MetaCompleteRoundFixedBellCorrIneqSym(OurEBI, iniP_proj, iniChannel, ins, outs);
outalpha = out{1};
outchannel = out{2};
outPproj = out{3};
outbell = out{4};

inputstate = final_state(ini_state(outalpha), outchannel);
val=evaluate_bell_ineq(EBIbellcoeffs, inputstate, outPproj, ins, outs);
prob = giveProbNDArray(inputstate, outPproj, ins, outs);
checkThatProbSumsToOne(prob,ins,outs);
disp([val,outalpha*4*sqrt(3)]);
% disp(ClassicalOptInequality2(outbell,ins,outs));