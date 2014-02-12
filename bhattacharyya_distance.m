function bdist = bhattacharyya_distance(histogram1, histogram2)
% This function computes the bhattacharyya distance between two normalised
% histograms.

bcoeff = sqrt(histogram1(:)')*sqrt(histogram2(:));
bdist = sqrt(1 - bcoeff);
end