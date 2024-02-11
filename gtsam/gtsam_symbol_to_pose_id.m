% function [prefix, frame_id] = gtsam_symbol_to_pose_id(symbol)
% Given a int64 gtsam symbol, return the corresponding tuple of prefix
% (char) and index (int)
function [robot_id, frame_id] = gtsam_symbol_to_pose_id(symbol)

symbol = int64(symbol);
gtsam_offsets = int64( ...
    [6989586621679009792; ...
    7061644215716937728; ...
    7133701809754865664; ...
    7205759403792793600; ...
    7277816997830721536; ...
    7349874591868649472; ...
    7421932185906577408; ...
    7493989779944505344; ...
    7566047373982433280; ...
    7638104968020361216]);
max_robot_id = size(gtsam_offsets, 1);
robot_id  = max_robot_id;

while robot_id > 0
    offset = gtsam_offsets(robot_id);
    if symbol >= offset
        frame_id = symbol - offset;
        return;
    else
        robot_id = robot_id - 1;
    end
end

error('Invalid input symbol: %d!', symbol);


end
