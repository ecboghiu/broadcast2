P = test1();
AA = randn(3)+sqrt(-1)*randn(3);
ia = find(triu(ones(3)));
P(AA(ia))

function [P]= test1()
X = sdpvar(3,3,'he','co')
A = sdpvar(3,3,'he','co')
Model = [X-A >= 0];
ia = find(triu(ones(3)));
P = optimizer(Model, trace(X),[],A(ia),X)
AA = randn(3)+sqrt(-1)*randn(3);
AA = AA + AA';
P(AA(ia))
optimize(X >= AA, trace(X))
value(X)
yalmip('clear');
end