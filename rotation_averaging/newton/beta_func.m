function beta = beta_func(theta, options)

if strcmp(options.rotation_distance, 'chordal')
    beta = sin(theta);
elseif strcmp(options.rotation_distance, 'geodesic')
    beta = theta / 2;
else
    error('Unknown rotation distance: %s', options.rotation_distance);
end

end