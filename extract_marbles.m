function [cc,cr,radii,n_detected, foreground]=extract_marbles(Imwork,Imback,fig1,fig2,fig3,fig15,index)
%Function extract_marbles
%Extract marbles from image
%marble_n is the n detected marble, and it is not yet a tracking id
%cc contains the centres of the marbles (<marble_n>,<frame>)
%cr contains the radii of the marbles (<marble_n>,<frame>)
%radii is an array containing the radius of each marble as is indexed by
%(<marble_n)
%n_detected is the number of detected marbles
% n_detected = 0 if failure, otherwise the number of detected marbles
  
  cc = 0;
  cr = 0;
  radii=0;
  n_detected=0;
  min_radius=3;
  max_radius=30;
  iThresholdBackground = 0.2;
  
  %%% Background subtractions for RGB, greyscale and Chromaticity:
  
  %% Greyscale background subtraction
  grey_back = sum(Imback,3)./3;
  grey_work = sum(Imwork,3)./3;
  greySubtraction = abs(grey_work - grey_back);
  
  % Adaptive threshold
  [counts, levels] = imhist(uint8(greySubtraction), 500);
  cdf = cumsum(counts)/sum(counts);
  GreyThreshold = max(10, levels( find(cdf < 0.97, 1, 'last') ));
  if isempty(GreyThreshold)
      GreyThreshold = 5;
  end
  grey_thresh = greySubtraction > GreyThreshold;
  cleaned = bwmorph(grey_thresh, 'clean', 1);
  cleaned = bwmorph(cleaned, 'close', 2);
  
  %% RGB background subtraction
  redSubtraction = abs(Imwork(:,:,1)-Imback(:,:,1));
  greenSubtraction = abs(Imwork(:,:,2)-Imback(:,:,2));
  blueSubtraction = abs(Imwork(:,:,3)-Imback(:,:,3));
  
  % Adaptive threshold for each colour channel
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
  
  %% Chromaticity backgorund subtraction
  Imback_ch = convertToChromaticity(Imback);
  Imwork_ch = convertToChromaticity(Imwork);
  
  RCHSubtraction = abs(Imwork_ch(:,:,1)-Imback_ch(:,:,1));
  GCHSubtraction = abs(Imwork_ch(:,:,2)-Imback_ch(:,:,2));
  BCHSubtraction = abs(Imwork_ch(:,:,3)-Imback_ch(:,:,3));
  
  % Adaptive threshold for each chromaticity channel
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
  
  redDetected = selectAreaRange(redDetected, min_radius, max_radius);
  greenDetected = selectAreaRange(greenDetected, min_radius, max_radius);
  blueDetected = selectAreaRange(blueDetected, min_radius, max_radius);
  RCHDetected = selectAreaRange(RCHDetected, min_radius, max_radius);
  GCHDetected = selectAreaRange(GCHDetected, min_radius, max_radius);
  BCHDetected = selectAreaRange(BCHDetected, min_radius, max_radius);
  greyDetected = selectAreaRange(greyDetected, min_radius, max_radius);
  
  % Combine all channels to find the final detections
  foreground = (redDetected | greenDetected | greenDetected ...
               | RCHDetected | GCHDetected | BCHDetected) ...
               & greyDetected;
  
  %{
  masked_rgb = (Imwork./255).*repmat(grey_thresh,1,1,3);
  masked_grey = sum(masked_rgb,3)./3;
  
  foreground = masked_grey;
  
  Imback_ch = convertToChromaticity(Imback);
  Imwork_ch = convertToChromaticity(Imwork);
  
  ch_subtract = abs( sum(Imback_ch(:,:,1:2),3)./2 - sum(Imwork_ch(:,:,1:2),3)./2 );
  %}

  foreground = bwmorph(foreground,'clean');
  foreground = bwmorph(foreground,'thicken');
  foreground = bwmorph(foreground,'close',10);
  foreground = bwmorph(foreground,'open',10);

  %Show the detections and the true image
  figure(1);
  imshow(foreground);
  figure(2);
  imshow(uint8(Imwork));
  
  

  % select labeled objects in a matrix
%  labeled = bwlabel(foremm,4);
  connected_components= bwconncomp(foreground,4);
%  stats = regionprops(labeled,['basic']);
%  [n_detected,W] = size(stats);

% make sure that there are marbles, which can have a radius between min
  % and max defined
  %put in im_detected_objects only objects matrix
  %code snippet adapted from Matlab13 help for regionprops function.
  stats =regionprops(connected_components,['basic']);
  idx_objects=find([stats.Area]>=(pi*(min_radius^2)) & [stats.Area]<=(pi*(max_radius^2)));
  im_detected_objects=ismember(labelmatrix(connected_components),idx_objects);
  stats_detected_objects=stats([stats.Area]>=(pi*(min_radius^2)) & [stats.Area]<=(pi*(max_radius^2)));

%  stats = stats(stats.Area>(pi*(min_radius^2)) & stats.Area<(pi*(max_radius^2)));
%  n_detected = [bwconncomp(im_detected_objects,4)].NumObjects;
  if isempty(idx_objects)
    n_detected=0;
    return   
  end
  
%{
  % do bubble sort (large to small) on regions in case there are more than 1

  id = zeros(N);
  for i = 1 : N
    id(i) = i;
  end
  for i = 1 : N-1
    for j = i+1 : N
      if stats(i).Area < stats(j).Area
        tmp = stats(i);
        stats(i) = stats(j);
        stats(j) = tmp;
        tmp = id(i);
        id(i) = id(j);
        id(j) = tmp;
      end
    end
  end
  %}


%  selected = (labeled==id(1));
  if fig3 > 0
    figure(fig3)
    clf
    imshow(im_detected_objects)
  end
%{
  % get center of mass and radius of largest
  centroid = stats(1).Centroid;
  radius = sqrt(stats(1).Area/pi);
  cc = centroid(1);
  cr = centroid(2);
  n_detected = 1;
  %}
  return
  
end

function im_detected_objects = selectAreaRange(bwimage, min_radius, max_radius)
  connected_components= bwconncomp(bwimage,4);
  stats = regionprops(connected_components,['basic']);
  idx_objects = find([stats.Area]>=(pi*(min_radius^2)) & [stats.Area]<=(pi*(max_radius^2)));
  im_detected_objects = ismember(labelmatrix(connected_components),idx_objects);
end
