% function R = euler_to_rotm2d(th)
% construct the 2D rotation matrix corresponding to input angle (radian)
function R = euler_to_rotm2d(th)

R = [cos(th), -sin(th); 
     sin(th), cos(th)];

end