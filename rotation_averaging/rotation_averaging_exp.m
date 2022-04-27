% R = rotation_averaging_exp(R, x, options)
% Given a current multiple rotation estimate R, apply the tangent vector x
% via exponential mapping to obtain the updated multiple rotation estimate.
% The input tangent vector x = vec([x1 ... xn]), where each xi is the
% rotation tangent vector xi.
%
% Yulun Tian
function R = rotation_averaging_exp(R, x, options)
    d = size(R,1);
    n = size(R,2) / d;
    if d == 2
        p = 1;  % intrisic dimension of a single pose
        Rexp = @exp2d;
    elseif d == 3
        p = 3; 
        Rexp = @exp3d;
    else
        error('Invalid dimension %i', d);
    end
    assert(size(x,1) == p*n);
    assert(size(x,2) == 1);
    
    for i = 1:n
        idxs = ((i-1)*d+1) : i*d;
        idxs_ = (i-1) *p +1: i * p;
        Ri = R(:, idxs);
        dRi = x(idxs_);
        if strcmp(options.tangent_space_parametrization, 'global')
            Ri_new = Rexp(dRi) * Ri;
        elseif strcmp(options.tangent_space_parametrization, 'local')
            Ri_new = Ri * Rexp(dRi);
        else
            error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization);
        end
        R(:, idxs) = Ri_new;
    end
end