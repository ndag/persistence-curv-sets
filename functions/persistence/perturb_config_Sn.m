function Y = perturb_config_Sn(X, sigma)
    % Perturb X by a normal epsilon with mean 0 and standard dev sigma
    Y = X + normrnd(0,sigma, size(X));
    
    % Renormalize Y so that it lands in Sn
    nor = sum(Y.^2,2).^0.5;                 % get norm
    Y = Y./(nor*ones(1, size(Y,2)));        % normalize
end