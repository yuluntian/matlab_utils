% function th = chordal2angular_3d(d)
% Given chordal distance, return the corresponding angular distance
function th = chordal2angular_3d(d)

th = 2 * asin(d/(2*sqrt(2)));

end