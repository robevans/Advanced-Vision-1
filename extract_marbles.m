function [centres,radii,n_detected]=extract_marbles(Imwork,Imback,fig1,fig2,fig3,fig15,index)
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
%       centres contains the centres of the marbles in x, y coordinates
%       radii is an array containing the radius of each marble as is indexed by
%       n_detected is the number of detected marbles
%       n_detected is 0 if failure, otherwise the number of detected marbles
  
  centres = 0;
  radii=0;
  n_detected=0;
  min_radius=5;
  max_radius=25;
  max_eccentricity=0.5;
  iThresholdBackground = 0.2;
  
  grey_back = sum(Imback,3)./3;
  grey_work = sum(Imwork,3)./3;
  
  grey_subtract = abs(grey_work - grey_back);
  
  % Adaptive threshold
  [counts, levels] = imhist(uint8(grey_subtract), 500);
  cdf = cumsum(counts)/sum(counts);
  threshold = max(10, levels( find(cdf < 0.98, 1, 'last') ))
  if isempty(threshold)
      threshold = 5;
  end
  grey_thresh = grey_subtract > threshold;
  
  cleaned = bwmorph(grey_thresh, 'clean', 1);
  cleaned = bwmorph(grey_thresh, 'close', 2);
  
  imshow(cleaned);
  
  masked_rgb = (Imwork./255).*repmat(grey_thresh,1,1,3);
  masked_grey = sum(masked_rgb,3)./3;
  
  foreground = masked_grey;
  
  Imback_ch = convertToChromaticity(Imback);
  Imwork_ch = convertToChromaticity(Imwork); 
  
  ch_subtract = abs( sum(Imback_ch(:,:,1:2),3)./2 - sum(Imwork_ch(:,:,1:2),3)./2 );
  
  [MR,MC,Dim] = size(Imback);

  % subtract background & select pixels with a big difference
 % fore = zeros(MR,MC);
  fore = (abs(Imwork(:,:,1)-Imback(:,:,1)) > iThresholdBackground) ...
     | (abs(Imwork(:,:,2) - Imback(:,:,2)) > iThresholdBackground); % ...
     %| (abs(Imwork(:,:,3) - Imback(:,:,3)) > iThresholdBackground);
  if fig15 > 0
    figure(fig1)
    clf
    %imshow(fore)
    %eval(['imwrite(uint8(fore),''BGONE/nobg',int2str(index),'.jpg'',''jpg'')']);  
  end

  % erode to remove small noise
  %foremm = bwmorph(fore,'erode',1);
  foremm = fore;
  foremm = bwmorph(foremm,'clean');
  foremm = bwmorph(foremm,'thicken');
  foremm = bwmorph(foremm,'close',10);
  foremm = bwmorph(foremm,'open',10);
  %foremm = bwmorph(foremm,'erode',1);
  %foremm = bwmorph(foremm,'thicken',5);
  %foremm = bwmorph(foremm,'erode',1);

  foreground=cleaned;
  if fig2 > 0
    %figure(fig2)
    %clf
    %imshow(foreground)
    imwrite(uint8(foreground),['BGCLEAN/',int2str(index),'.jpg'],'jpg');  
  end
  


  % select labeled objects in a matrix
%  labeled = bwlabel(foremm,4);
  connected_components= bwconncomp(foreground,4);
%  stats = regionprops(labeled,['basic']);
%  [n_detected,W] = size(stats);

% make sure that there are marbles, which can have a radius between min
  % and max defined
  %put in im_detected_objects only objects matrix
  %code snippet adapted from Matlab13 help for regionprops function.
  stats =regionprops(connected_components,'Eccentricity','Centroid','Area');
  idx_objects=find([stats.Area]>=(pi*(min_radius^2)) & [stats.Area]<=(pi*(max_radius^2)) & [stats.Eccentricity] < max_eccentricity);  
  im_detected_objects=ismember(labelmatrix(connected_components),idx_objects);
  stats_detected_objects=stats([stats.Area]>=(pi*(min_radius^2)) & [stats.Area]<=(pi*(max_radius^2)) & [stats.Eccentricity] < max_eccentricity );

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
%    imwrite(uint8(im_detected_objects),['BGDETECT/',int2str(index),'.jpg'],'jpg');      
  end

  centres=zeros(2,length(stats_detected_objects));
  radii=zeros(length(stats_detected_objects),1);
  % get center of mass and radius of largest
  for i=1:length(stats_detected_objects)
      centres(:,i)=stats_detected_objects(i).Centroid;
      radii(i)=sqrt(stats_detected_objects(i).Area/pi);
      n_detected=n_detected+1;
  end
  
  return
