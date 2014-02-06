
function [ Images_marbles,success ] = detect_marbles( directory )
%Function detect_marbles
%AV1 Practical 1
% This function detects marbles in images located in the provided directory
% and outputs its detection in an array, it returns too the number of
% marbles detected
% Parts of code based detectp.m
% Input :  directory - The directory where images files are located in
% format %d.jpg
% Output : Images_marbles - An array of images detected
%          sucess - Number of detected marbles, 0 if none which corresponds
%          to failure, non-zero is sucess
%Robert Evans, Francisco Aguirre
%Temporary
Images_marbles=true;
success=true;

imgBackGround = double(imread(strcat(directory,'1.','jpg')));
[MR,MC,Dim] = size(imgBackGround);

%Counting how many data images we got
num_Images=length(dir('SEQ1/*jpg'));

 % loop over all images
fig1=0;
fig2=0; %show image after erosion
fig15=0; %debug, show foreground minus background
fig3=3; %debug, show detected objects
fig4=0;

%num_Images=10;

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
  
  [centres(:,:,i),radii(:,i),n_detected(i)]=extract_marbles(Imwork,imgBackGround,fig1,fig2,fig3);

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

