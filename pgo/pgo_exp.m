% function [R, t] = pgo_exp(R, t, x, options)
% Given a current pose graph estimate (R, t), apply the tangent vector x
% via exponential mapping to obtain the updated pose graph estimate.
% The input tangent vector x = vec([x1 ... xn]), where each xi is the pose
% tangent vector xi = [xi_rotation xi_translation]
%
% Yulun Tian
function [R, t] = pgo_exp(R, t, x, options)
    d = size(t,1);
    n = size(t,2);
    if d == 2
        p = 3;  % intrisic dimension of a single pose
        Rexp = @exp2d;
    elseif d == 3
        p = 6; 
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
        ti = t(:,i);
        xi = x(idxs_);
        if d == 2
            dRi = xi(1);
            dti = xi(2:3);
        else
            dRi = xi(1:3);
            dti = xi(4:6);
        end
        ti_new = ti + dti;
        if strcmp(options.tangent_space_parametrization, 'global')
            Ri_new = Rexp(dRi) * Ri;
        elseif strcmp(options.tangent_space_parametrization, 'local')
            Ri_new = Ri * Rexp(dRi);
        else
            error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization);
        end
        R(:, idxs) = Ri_new;
        t(:, i) = ti_new;
    end

end