x = sdpvar(1);
tic
optimize(x>= 0, x,sdpsettings('debug',1,'verbose',0));
toc