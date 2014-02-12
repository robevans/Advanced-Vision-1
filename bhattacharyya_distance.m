function bdist = bhattacharyya_distance(histogram1, histogram2)
% This function computes the bhattacharyya distance between two normalised
% histograms.

% Computes the bhattacharyya coefficient
bcoeff = sqrt(histogram1(:)')*sqrt(histogram2(:));

% Uses the coefficient to compute the bhattacharyya distance.
bdist = sqrt(1 - bcoeff);
end