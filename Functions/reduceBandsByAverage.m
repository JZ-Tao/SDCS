function reducedImg = reduceBandsByAverage(x, targetBands)
    % x: Input multi-band image
    % targetBands: Number of target bands, e.g. 4

    [height, width, bands] = size(x); % Get the size and number of bands of the original image
    bandGroups = bands / targetBands; % Calculate the number of original bands to be merged for each target band

    % Pre-allocated memory space for storing the final image
    reducedImg = zeros(height, width, targetBands, class(x));

    for i = 1:targetBands
        % Calculate the original band range corresponding to the current target band
        startBand = round((i-1) * bandGroups) + 1;
        endBand = round(i * bandGroups);

        % Calculate the average of these bands
        reducedImg(:,:,i) = mean(x(:,:,startBand:endBand), 3);
    end
end