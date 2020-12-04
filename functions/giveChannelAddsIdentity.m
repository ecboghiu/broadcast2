function output_choi_channel = giveChannelAddsIdentity(dimInput,dimIdentity, placement)
    % 'placement' should be either "left" or "right"
    output_choi_channel = 0;
    for i=1:dimInput
        for j=1:dimInput
            ketbra_ij = ketbra(i,j,dimInput);
            if placement == "left"
                output_choi_channel = output_choi_channel + kron(ketbra_ij, kron(eye(dimIdentity)/dimIdentity, ketbra_ij));
            elseif placement == "right"
                output_choi_channel = output_choi_channel + kron(ketbra_ij, kron(ketbra_ij, eye(dimIdentity)/dimIdentity));
            else
                error('incorrect "placement" value: should be either "left" or "right"');
            end
        end
    end
end

