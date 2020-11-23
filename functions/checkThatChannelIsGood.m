function result = checkThatChannelIsGood(channel, dimx, dimy)

choi = ChoiMatrix(channel);

tracedoutput = PartialTrace(choi, [2], [dimx, dimy]);
diff = tracedoutput - eye(dimx);
assert( norm(diff) <= 1e-6, 'Choi matrix condition not satisfied.');

choieigs = eig(choi);
assert( all(choieigs >= -1e-6), 'Choi matrix not semi-definite positive.');


result=true;
end