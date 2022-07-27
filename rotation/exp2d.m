% function R = exp2d(w)
% Compute the exponential mapping of SO(2). Input w is a scalar that
% corresponds to an element of Lie algebra so(2)
% 
% Yulun Tian
function R = exp2d(w)

% W = w * gen2;
% assert(size(W,1) == 2);
% assert(size(W,2) == 2);
% R = expm(W);
R = [cos(w) -sin(w);
       sin(w)   cos(w)];

end