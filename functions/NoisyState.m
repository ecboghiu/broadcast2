function [state] = NoisyState(noise,params)
state_family_name = params.name;

if strcmp(state_family_name, 'partially_entangled')
    % partially entangled
    chi = params.CONST_CHI;
    id_on = params.IDENTITY_PLACEMENT;
    state = NoisyPartiallyEntangled(noise, chi, id_on);
elseif strcmp(state_family_name, 'werner')
    state = NoisyWernerState(noise);
else
   error("Invalid input name"); 
end

end

