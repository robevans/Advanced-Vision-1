function imageChromaticity = convertToChromaticity(imageRGB)
% Converts to chromaticity
imageChromaticity = zeros(size(imageRGB));
rgbsum = sum(imageRGB,3);
imageChromaticity(:,:,1) = double(imageRGB(:,:,1)) ./ rgbsum;
imageChromaticity(:,:,2) = double(imageRGB(:,:,2)) ./ rgbsum;
imageChromaticity(:,:,3) = double(imageRGB(:,:,3)) ./ (rgbsum/3);
imageChromaticity(:,:,4) = double(sum(imageRGB(:,:,1:2),3)) ./ rgbsum;