%showgt(1,2,71)
% 
% This function plots the centres of the marbles in 2 frames
% allfignum: shows the connected trajectories of all marbles in the same frame
% tmpfignum: shows the detected marble centre for each frame
% numframes: number of video frames in the sequence
%
% You have to edit the lines for variables imagefile_name and gtfile_name
%    for the specific file and path names
%
% Type return after each frame to afvance to the next
function showgt(allfignum,tmpfignum,numframes)
    
%imagefile_name='/afs/inf.ed.ac.uk/user/r/rbf/dept/ADVANCED_VISION/PRACTICALS/PRACT1314/PRACT1/SEQ1/'; % full path to groundtruth
imagefile_name='./SEQ1/'; % relative path to groundtruth and images

gtfile_name = 'gtSeq1.mat';
file_format='.jpg';
load(gtfile_name) % load ground truth

% plot all marbles found in current frame
boxheight=1; % half of drawn box height
boxwidth=1;
figure(tmpfignum)   % figure to overlap current detection
for frame = 1 : numframes  % loop over all frames
  current_frame=imread([imagefile_name num2str(frame) file_format]);
  imshow(current_frame);
  hold on;

  for marblenum=1:size(new_marbles_comingFromRight,2)       % loop over the individual marbles coming from the right
    index_a=find(new_marbles_comingFromRight(marblenum).frame_numbers(:)==frame); % get index in full framelist for current frame
    if ~isempty(index_a)
      x=new_marbles_comingFromRight(marblenum).row_of_centers(index_a); % x left most pixel
      y=new_marbles_comingFromRight(marblenum).col_of_centers(index_a); % y of left most pixel
[1,frame,marblenum,x,y]
      plot([x-boxwidth,x-boxwidth,x+boxwidth,x+boxwidth,x-boxwidth],[y-boxheight,y+boxheight,y+boxheight,y-boxheight,y-boxheight],'r-') % draw box around centre
    end
  end
  for marblenum=1:size(new_marbles_comingFromLeft,2)       % loop over the individual marbles coming from the left
    index_a=find(new_marbles_comingFromLeft(marblenum).frame_numbers(:)==frame); % get index in full framelist for current frame
    if ~isempty(index_a)
      x=new_marbles_comingFromLeft(marblenum).row_of_centers(index_a); % x left most pixel
      y=new_marbles_comingFromLeft(marblenum).col_of_centers(index_a); % y of left most pixel
[2,frame,marblenum,x,y]
      plot([x-boxwidth,x-boxwidth,x+boxwidth,x+boxwidth,x-boxwidth],[y-boxheight,y+boxheight,y+boxheight,y-boxheight,y-boxheight],'r-') % draw box around centre
    end
  end
  drawnow;
end

pause


% plot all marble trajectories on frame 1
figure(allfignum)   % figure to show all trajectories
current_frame=imread([imagefile_name '1' file_format]);
imshow(current_frame);
hold on;

numdetects=zeros(9,1);
allx=zeros(9,100);
ally=zeros(9,100);
colourlist=['wrgbykmc'];

for marblenum=1:size(new_marbles_comingFromRight,2)       % loop over the individual marbles coming from the right
    frameList=new_marbles_comingFromRight(marblenum).frame_numbers(:);  % set of frame this marble appears in
    if ~isempty(frameList)
      last_frame=max(frameList);
      first_frame=min(frameList);
      for frame=first_frame:last_frame   % loop over all frames in list
        index_a=find(frameList(:)==frame); % get index in full framelist for current frame
        numdetects(marblenum)=numdetects(marblenum)+1;
        allx(numdetects(marblenum))=new_marbles_comingFromRight(marblenum).row_of_centers(index_a); % x left most pixel
        ally(numdetects(marblenum))=new_marbles_comingFromRight(marblenum).col_of_centers(index_a); % y of left most pixel
      end
[1,marblenum,numdetects(marblenum)]
allx(1:numdetects(marblenum))
      plot(allx(1:numdetects(marblenum)),ally(1:numdetects(marblenum)),[colourlist(mod(marblenum,8)+1) '*-'])
    end
end
pause

numdetects=zeros(9,1);
for marblenum=1:size(new_marbles_comingFromLeft,2)       % loop over the individual marbles coming from the left
    frameList=new_marbles_comingFromLeft(marblenum).frame_numbers(:);  % set of frame this marble appears in
    if ~isempty(frameList)
      last_frame=max(frameList);
      first_frame=min(frameList);
      for frame=first_frame:last_frame   % loop over all frames in list
        index_a=find(frameList(:)==frame); % get index in full framelist for current frame
        numdetects(marblenum)=numdetects(marblenum)+1;
        allx(numdetects(marblenum))=new_marbles_comingFromLeft(marblenum).row_of_centers(index_a); % x left most pixel
        ally(numdetects(marblenum))=new_marbles_comingFromLeft(marblenum).col_of_centers(index_a); % y of left most pixel
      end
[2,marblenum,numdetects(marblenum)]
allx(1:numdetects(marblenum))
      plot(allx(1:numdetects(marblenum)),ally(1:numdetects(marblenum)),[colourlist(mod(marblenum,8)+1) '*-'])
    end
end
