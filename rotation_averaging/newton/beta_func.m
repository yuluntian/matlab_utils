function beta = beta_func(theta, options)

if strcmp(options.rotation_distance, 'chordal')
    beta = sin(theta);
else
    error('Unknown rotation distance: %s', options.rotation_distance);
end

end