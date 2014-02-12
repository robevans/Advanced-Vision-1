function bdist = bhattacharyya_distance(histogram1, histogram2)
% This function computes the bhattacharyya distance between two normalised
% histograms.

%bcoeff = sqrt(histogram1(:)')*sqrt(histogram2(:));

% get number of bins
% (i.e. dimension 2 from Nx1 inputs)
bins = size(histogram1, 2);

% estimate the bhattacharyya co-efficient
bcoeff = 0;

for i=1:bins
    bcoeff = bcoeff + sqrt(histogram1(i) * histogram2(i));    
end

% get the distance between the two distributions as follows
bdist = sqrt(1 - bcoeff);
end