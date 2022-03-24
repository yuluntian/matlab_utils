function gamma = gamma_func(theta, options)

if strcmp(options.rotation_distance, 'chordal')
    if abs(theta) < 1e-12
        gamma = 0;
    else
        gamma = 2*cos(theta) - alpha_func(theta, options);
    end
else
    error('Unknown rotation distance: %s', options.rotation_distance);
end

end