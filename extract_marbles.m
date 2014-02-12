function [matMarbles,n_detected,new_background]=extract_marbles(Imwork,Imback,fig1,fig2,fig3,bShowImages)
%Function extract_marbles
%AV Practical 1
%Extract marbles from image
%Input: Imwork - Current image to extract from
%       Imback - Background image
%       fig1   - if set, shows background substraction, usually set to
%       value 1
%       fig2   - if set, shows foreground, usually set to value 2
%       fig3   - if set, shows detected objects, usually set to value 3
%       fig15  - used for debugging, usually set to value 15. Not currently
%       used
%       index  - Current image index
%Output:
%       matCentres contains the centres of the marbles in x, y coordinates
%       plus the radii
%       radii is an array containing the radius of each marble as is indexed by
%       n_detected is the number of detected marbles
%       n_detected is 0 if failure, otherwise the number of detected marbles

%matMarbles contains x, y, radius and the histogram of the marble, 256
%values for the later.
matMarbles = zeros(1,3+256);
n_detected=0;
min_radius=3;
max_radius=30;
max_eccentricity=0.9;
max_circularity=2;

%%% Independent Background subtractions and adaptive thresholding
%   for RGB, greyscale and Chromaticity channels:

%% Greyscale background subtraction
grey_back = sum(Imback,3)./3;
grey_work = sum(Imwork,3)./3;
greySubtraction = abs(grey_work - grey_back);

%% RGB background subtraction
redSubtraction = abs(Imwork(:,:,1)-Imback(:,:,1));
greenSubtraction = abs(Imwork(:,:,2)-Imback(:,:,2));
blueSubtraction = abs(Imwork(:,:,3)-Imback(:,:,3));

%% Chromaticity backgorund subtraction
Imback_ch = convertToChromaticity(Imback);
Imwork_ch = convertToChromaticity(Imwork);

RCHSubtraction = abs(Imwork_ch(:,:,1)-Imback_ch(:,:,1));
GCHSubtraction = abs(Imwork_ch(:,:,2)-Imback_ch(:,:,2));
BCHSubtraction = abs(Imwork_ch(:,:,3)-Imback_ch(:,:,3));

%% Adaptive threshold for greyscale subtraction
[counts, levels] = imhist(uint8(greySubtraction), 500);
cdf = cumsum(counts)/sum(counts);
GreyThreshold = max(10, levels( find(cdf < 0.97, 1, 'last') ));
if isempty(GreyThreshold)
    GreyThreshold = 5;
end

%% Adaptive threshold for RGB
[Rcounts, Rlevels] = imhist(uint8(redSubtraction), 500);
[Gcounts, Glevels] = imhist(uint8(greenSubtraction), 500);
[Bcounts, Blevels] = imhist(uint8(blueSubtraction), 500);
Rthreshold = max(10, Rlevels( find(cumsum(Rcounts)/sum(Rcounts) < 0.99, 1, 'last') ));
Gthreshold = max(10, Glevels( find(cumsum(Gcounts)/sum(Gcounts) < 0.99, 1, 'last') ));
Bthreshold = max(1, Blevels( find(cumsum(Bcounts)/sum(Bcounts) < 0.99, 1, 'last') ));
if isempty(Rthreshold)
    Rthreshold = 5;
end
if isempty(Gthreshold)
    Gthreshold = 5;
end
if isempty(Bthreshold)
    Bthreshold = 5;
end

%% Adaptive threshold for each chromaticity channel
[RCHcounts, RCHlevels] = imhist(RCHSubtraction, 500);
[GCHcounts, GCHlevels] = imhist(GCHSubtraction, 500);
[satcounts, satlevels] = imhist(BCHSubtraction, 500);
RCHthreshold = max(0.01, RCHlevels( find(cumsum(RCHcounts)/sum(RCHcounts) < 0.99, 1, 'last') ));
GCHthreshold = max(0.01, GCHlevels( find(cumsum(GCHcounts)/sum(GCHcounts) < 0.99, 1, 'last') ));
BCHthreshold = max(0.01, satlevels( find(cumsum(satcounts)/sum(satcounts) < 0.99, 1, 'last') ));
if isempty(RCHthreshold)
    RCHthreshold = 0.01;
end
if isempty(GCHthreshold)
    GCHthreshold = 0.01;
end
if isempty(BCHthreshold)
    BCHthreshold = 0.01;
end

%% Process the detections for each channel
redDetected =redSubtraction > Rthreshold;
greenDetected = greenSubtraction > Gthreshold;
blueDetected = blueSubtraction > Bthreshold;
RCHDetected = RCHSubtraction > RCHthreshold;
GCHDetected = GCHSubtraction > GCHthreshold;
BCHDetected = BCHSubtraction > BCHthreshold;
greyDetected= greySubtraction > GreyThreshold;

% Clean up each channel
redDetected = bwmorph(redDetected, 'clean', 2);
greenDetected = bwmorph(greenDetected, 'clean', 2);
blueDetected = bwmorph(blueDetected, 'clean', 2);
RCHDetected = bwmorph(RCHDetected, 'clean', 2);
GCHDetected = bwmorph(GCHDetected, 'clean', 2);
BCHDetected = bwmorph(BCHDetected, 'clean', 2);
greyDetected= bwmorph(greyDetected, 'clean', 2);

redDetected = bwmorph(redDetected, 'close', 5);
greenDetected = bwmorph(greenDetected, 'close', 5);
blueDetected = bwmorph(blueDetected, 'close', 5);
RCHDetected = bwmorph(RCHDetected, 'close', 5);
GCHDetected = bwmorph(GCHDetected, 'close', 5);
BCHDetected = bwmorph(BCHDetected, 'close', 5);
greyDetected= bwmorph(greyDetected, 'close', 5);

% Reject objects too large or too small
redDetected = selectAreaRange(redDetected, min_radius, max_radius);
greenDetected = selectAreaRange(greenDetected, min_radius, max_radius);
blueDetected = selectAreaRange(blueDetected, min_radius, max_radius);
RCHDetected = selectAreaRange(RCHDetected, min_radius, max_radius);
GCHDetected = selectAreaRange(GCHDetected, min_radius, max_radius);
BCHDetected = selectAreaRange(BCHDetected, min_radius, max_radius);
greyDetected = selectAreaRange(greyDetected, min_radius, max_radius);

%% Combine detections from all channels.
foreground = (redDetected | greenDetected | blueDetected ...
    | RCHDetected | GCHDetected | BCHDetected) ...
    & greyDetected;

% Do some clean up on the composite detections.
foreground = bwmorph(foreground,'clean');
%foreground = bwmorph(foreground,'thicken');
foreground = bwmorph(foreground,'close',10);
foreground = bwmorph(foreground,'open',20);

%Show the detections and the true image
if fig1 > 0 && fig2 > 0 && bShowImages
    figure(fig1);
    imshow(foreground);
    figure(fig2);
    imshow(uint8(Imwork));
end

% Extract the detected background so it can be averaged into the current background.
new_background = Imwork.*repmat(~foreground,1,1,3);
if bShowImages
    figure(17);
    imshow(uint8(new_background));
end

% select labeled objects in a matrix
connected_components= bwconncomp(foreground,4);

% Filter detections by Area, eccentricity and circularity properties.
stats = regionprops(connected_components,'Eccentricity','Centroid','Area','Perimeter');
idx_objects = find([stats.Area]>=(pi*(min_radius^2)) & [stats.Area]<=(pi*(max_radius^2)) & [stats.Eccentricity] < max_eccentricity & (([stats.Perimeter] .^ 2) ./ (4 * pi * [stats.Area])) < max_circularity);

% Create final detections image
im_detected_objects = ismember(labelmatrix(connected_components),idx_objects);
stats_detected_objects = stats([stats.Area]>=(pi*(min_radius^2)) & [stats.Area]<=(pi*(max_radius^2)) & [stats.Eccentricity] < max_eccentricity );

% Return if nothing was detected.
if isempty(idx_objects)
    n_detected=0;
    return
end

% Display detections image.
if fig3 > 0 && bShowImages
    figure(fig3)
    clf
    imshow(im_detected_objects)
end

% Return detected objects' properties
matMarbles=[];
% get center of mass and radius of each object
for i=1:length(stats_detected_objects)
    matMarbles(i,1:2)=stats_detected_objects(i).Centroid;
    matMarbles(i,3)=sqrt(stats_detected_objects(i).Area/pi);
    %Get the histogram
    matMarbles(i,4:(256+3))=(histogramOfCircleAroundPoint(matMarbles(i,1),matMarbles(i,2),matMarbles(i,3),Imwork))';
    n_detected=n_detected+1;
end

return

end

% Helper function to control for area in individual channels.
function im_detected_objects = selectAreaRange(bwimage, min_radius, max_radius)
connected_components= bwconncomp(bwimage,4);
stats = regionprops(connected_components,'basic');
idx_objects = find([stats.Area]>=(pi*(min_radius^2)) & [stats.Area]<=(pi*(max_radius^2)));
im_detected_objects = ismember(labelmatrix(connected_components),idx_objects);
end
