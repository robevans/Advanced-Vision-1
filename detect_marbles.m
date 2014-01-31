%Function detect_marbles
%AV1 Practical 1
%Robert Evans, Francisco Aguirre

function [ Images_marbles,success ] = detect_marbles( directory )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Temporary
Images_marbles=true;

success=true;

%Code borrowed from detect.m
% load the background image
Imback = double(imread(strcat(directory,'1.','jpg')));
[MR,MC,Dim] = size(Imback);

%Counting how many data images we got
num_Images=length(dir('SEQ1/*jpg'));

% loop over all images
fig1=1;
fig2=0;
fig15=0;
fig3=2; %debug
fig4=4;
for i = 1 : num_Images
  % load image
  Im = imread([strcat(directory,int2str(i)), '.jpg'],'jpg');
  if fig1 > 0
    figure(fig1)
    clf
    imshow(Im)
  end
  Imwork = double(Im);

  %Extract marbles from image
  %marble_n is the n detected marble, and it is not yet a tracking id
  %cc contains the centres of the marbles (<marble_n>,<frame>)
  %cr contains the radii of the marbles (<marble_n>,<frame>)
  %radii is an array containing the radius of each marble as is indexed by
  %(<marble_n)
  %n_detected is the number of detected marbles
  
  [cc(:,i),cr(:,i),radii,n_detected]=extract_marbles(Imwork,Imback,fig1,fig2,fig3,fig15,i);
  if flag==0
    continue
  end
%{
  if fig1 > 0
    figure(fig1)
    hold on
    for c = -0.97*radius: radius/20 : 0.97*radius
      r = sqrt(radius^2-c^2);
      plot(cc(i)+c,cr(i)+r,'g.')
      plot(cc(i)+c,cr(i)-r,'g.')
    end
    %eval(['saveas(gcf,''TRACK/trk',int2str(i-1),'.jpg'',''jpg'')']);  
  end
%}
      pause(0.3)
end

% show positions
if fig4 > 0
  figure(fig4)
  hold on
  clf
  plot(cc,'r*')
  plot(cr,'g*')
end




end

