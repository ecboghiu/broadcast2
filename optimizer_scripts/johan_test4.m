R = sdpvar(3,3);
I = sdpvar(3,3,'skew');
Z = [R I;-I R];
obj = trace(R);
constr = [Z <= 1];
opt=optimize(constr,-obj,sdpsettings('solver','mosek'));
value(obj)

opt=optimizer(constr,-obj,[],[],R);
simple_numeric_skew = [[0,1,1];[-1,0,1];[-1,-1,0]];
R2=opt()
trace(R2)
