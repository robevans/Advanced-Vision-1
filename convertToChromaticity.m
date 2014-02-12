function imageChromaticity = convertToChromaticity(imageRGB)
% This function converts an RGB image to chromaticity coordinates, plus a
% fourth channel containing saturation.
imageChromaticity = zeros(size(imageRGB));
rgbsum = sum(imageRGB,3);
imageChromaticity(:,:,1) = double(imageRGB(:,:,1)) ./ rgbsum;
imageChromaticity(:,:,2) = double(imageRGB(:,:,2)) ./ rgbsum;
imageChromaticity(:,:,3) = double(imageRGB(:,:,3)) ./ rgbsum;
imageChromaticity(:,:,4) = rgbsum ./ 3;