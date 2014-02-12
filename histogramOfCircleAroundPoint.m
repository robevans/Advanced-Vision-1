function circle_histogram = histogramOfCircleAroundPoint(x,y,radius,RGBimage)
% This function extracts the pixel values around the given point and
% returns a normalised histogram of those values.

%Convert to chromaticity
chromaticity_image = convertToChromaticity(RGBimage);
RG_chromaticity_sum = sum(chromaticity_image(:,:,1:2),3) / 2;

% Extract circle pixels
[rows, cols, ~] = size(RG_chromaticity_sum);
[xx,yy] = ndgrid( (1:rows)-x, (1:cols)-y );
mask = (xx.^2 + yy.^2)<radius^2;

% Get histogram of circle pixels and normalise it
[h, c] = imhist(RG_chromaticity_sum(mask));
circle_histogram = normaliseHistogram(h, c);
end

function normalised_histogram = normaliseHistogram(histogram, binCenters)
% This function normalises a histogram to sum to 1.
dx = diff(binCenters(1:2));
normalised_histogram = histogram/sum(histogram*dx);
end