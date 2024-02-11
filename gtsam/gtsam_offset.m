% Compute the index offset for a robot in gtsam
function offset = gtsam_offset(robot_id)

if robot_id == 1
    % Add prefix 'a'
    offset = int64(6989586621679009792);
elseif robot_id == 2
    % Add prefix 'b'
    offset = int64(7061644215716937728);
elseif robot_id == 3
    % Add prefix 'c'
    offset = int64(7133701809754865664);
elseif robot_id == 4
    % Add prefix 'd'
    offset = int64(7205759403792793600);
elseif robot_id == 5
    % Add prefix 'e'
    offset = int64(7277816997830721536);
elseif robot_id == 6
    % Add prefix 'f'
    offset = int64(7349874591868649472);
elseif robot_id == 7
    % Add prefix 'g'
    offset = int64(7421932185906577408);
elseif robot_id == 8
    % Add prefix 'h'
    offset = int64(7493989779944505344);
elseif robot_id == 9
    % Add prefix 'i'
    offset = int64(7566047373982433280);
elseif robot_id == 10
    % Add prefix 'j'
    offset = int64(7638104968020361216);
else
    error('Unsupported robot id for gtsam: %i', robot_id);
end

end