function [cc,cr,radii,n_detected]=extract_marbles(Imwork,Imback,fig1,fig2,fig3,fig15,index)
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
  min_radius=10;
  max_radius=25;
  
  
  [MR,MC,Dim] = size(Imback);

  % subtract background & select pixels with a big difference
  fore = zeros(MR,MC);
  fore = (abs(Imwork(:,:,1)-Imback(:,:,1)) > 10) ...
     | (abs(Imwork(:,:,2) - Imback(:,:,2)) > 10) ...
     | (abs(Imwork(:,:,3) - Imback(:,:,3)) > 10);
  if fig15 > 0
    figure(fig15)
    clf
    imshow(fore)
    %eval(['imwrite(uint8(fore),''BGONE/nobg',int2str(index),'.jpg'',''jpg'')']);  
  end

  % erode to remove small noise
  foremm = bwmorph(fore,'erode',2);

  if fig2 > 0
    figure(fig2)
    clf
    imshow(foremm)
    %eval(['imwrite(uint8(foremm),''BCLEAN/cln',int2str(index),'.jpg'',''jpg'')']);  
  end

  % select labeled objects in a matrix
%  labeled = bwlabel(foremm,4);
  connected_components= bwconncomp(foremm,4);
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
