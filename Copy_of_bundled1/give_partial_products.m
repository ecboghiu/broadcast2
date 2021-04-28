function [partial_products] = give_partial_products(povms, bellcoeffs, ins, outs)

nrparties = length(ins);
partial_products=cell(3,1);

for party_to_ignore=1:nrparties
    partial_products{party_to_ignore} = cell(ins(party_to_ignore), outs(party_to_ignore));
    allbutone = logical(ones(nrparties,1)); 
    allbutone(party_to_ignore) = false;
    ins_without_party_to_ignore = ins(allbutone);
    outs_without_party_to_ignore = outs(allbutone);
    
    dimen = prod(outs_without_party_to_ignore);
    for x=1:ins(party_to_ignore)
        for a=1:outs(party_to_ignore)
            partial_products{party_to_ignore}{x,a} = 0;
            for y=1:ins_without_party_to_ignore(1)
                for z=1:ins_without_party_to_ignore(2)
                    for b=1:outs_without_party_to_ignore(1)
                        for c=1:outs_without_party_to_ignore(2)
                            if party_to_ignore == 1
                                term = bellcoeffs(x,y,z,a,b,c) * kron(povms{2}{y}{b},povms{3}{z}{c});
                            elseif party_to_ignore == 2
                                term = bellcoeffs(y,x,z,b,a,c) * kron(povms{1}{y}{b},povms{3}{z}{c});
                            elseif party_to_ignore == 3
                                term = bellcoeffs(y,z,x,b,c,a) * kron(povms{1}{y}{b},povms{2}{z}{c});
                            else
                               error('more than 3 parties not supported'); 
                            end
                            partial_products{party_to_ignore}{x,a} = partial_products{party_to_ignore}{x,a} + term;
                        end
                    end
                end
            end
        end
    end
    
end

end