function alpha = alpha_func(theta, options)

if strcmp(options.rotation_distance, 'chordal')
    if abs(theta) < 1e-12
        alpha = 2;
    else
        alpha = sin(theta) * cot(theta/2);
    end
elseif strcmp(options.rotation_distance, 'geodesic')
    if abs(theta) < 1e-12
        alpha = 1;
    else
        alpha = 0.5 * theta * cot(theta/2);
    end
else
    error('Unknown rotation distance: %s', options.rotation_distance);
end

end