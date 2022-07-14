% Implemented according to SE-Sync Appendix A: 
% https://arxiv.org/pdf/1612.07386.pdf
function R = Langevin_sampler_2D(M, kappa)

% sample a rotation angle in vM(0, 2*kappa)
theta = circ_vmrnd(0, 2*kappa, 1);

% construct corresponding rotation matrix as in eq. (54)
P = [cos(theta) -sin(theta); ...
     sin(theta) cos(theta)];

R = M*P;

end
