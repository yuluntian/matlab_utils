% function Gi = gen3(i)
% Return the generator matrix on so(3)
% 
% Yulun Tian
function Gi = gen3(i)

if i == 1
    Gi = [0 0 0; 0 0 -1; 0 1 0];
elseif i == 2
    Gi = [0 0 1; 0 0 0; -1 0 0];
elseif i == 3
    Gi = [0 -1 0; 1 0 0; 0 0 0];
else
    error('Invalid input: %i!',i);
end

end