% function d = angular2chordal(th)
% Convert angular distance in radian to chordal distance
function d = angular2chordal(th)

d = 2*sqrt(2)*sin(th/2);

end