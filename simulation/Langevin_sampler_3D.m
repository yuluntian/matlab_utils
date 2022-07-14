% Implemented according to SE-Sync Appendix A: 
% https://arxiv.org/pdf/1612.07386.pdf
function R = Langevin_sampler_3D(M, kappa)

% sample a rotation angle in vM(0, 2*kappa)
theta = circ_vmrnd(0, 2*kappa, 1);

% sample an uniformly distributed axis
S2_manifold = spherefactory(3);
vhat = S2_manifold.rand();
vhat = vhat / norm(vhat);
vhat_mat = [0 -vhat(3) vhat(2);
             vhat(3) 0 -vhat(1);
             -vhat(2) vhat(1) 0];

P = expm(theta*vhat_mat);
R = M*P;


end
