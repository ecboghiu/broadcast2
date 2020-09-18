function corr = ToCorrelatorNotationWithSymProbInput(syminput)

nameinput = char(syminput);

% assuming input of the form p_abc_xyz and no integer greater than 10

a = str2double(nameinput(3));
b = str2double(nameinput(4));
c = str2double(nameinput(5));
x = str2double(nameinput(7));
y = str2double(nameinput(8));
z = str2double(nameinput(9));


corr = ToCorrelatorNotationWithInOut([x,y,z],[a,b,c]);
end