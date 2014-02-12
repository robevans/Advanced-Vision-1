function [ success ] = ground_truth_statistics( directory, Tracking_marbles, marblesWithIDs )
%Function ground_truth_statistics
%AV 1 Practical 1 20140204
%This function compares ground truth statistics in getSeq1.mat file
%to statistics obtained by the program as output from tracking marbles.
%
%Input: directory - Directory where ground file is located
%
%       Tracking_marbles - A struct containing frames and for each frame
%       and (x,y) coordinates of detected marbles
%
%       marblesWithIDs - A matrix with centre x,y and marble ids for each
%       frame.
%
%Output: success - True if there was no error
%
%Portions of code come from showgt where the ground truth file is used to
%display trajectories

%Ground truth file
gtfile_name = 'gtSeq1.mat';

%radius for detection
iDetectionRadius=10;

%Counting how many data images we got
num_Images=length(dir(strcat(directory,'*jpg')));

%Matrix to hold statistics, per frame, marbles in ground truth and marbles detected
%within 10 pixels
matGroundTruthStats=zeros(num_Images,2);

sum_of_distances_from_ground_truth = 0;

% load ground truth
%This file loads the following variables
%new_marbles_comingFromRight
%new_marbles_comingFromLeft
if ~exist(gtfile_name,'file')
    success=0;
    return;
end

load(gtfile_name);


for frame = 1 : num_Images  % loop over all frames
    
    %% Count number of detections within 10 pixels of ground truth, and average distance of those detections.
    for marblenum=1:size(new_marbles_comingFromRight,2)       % loop over the individual marbles coming from the right
        index_a=find(new_marbles_comingFromRight(marblenum).frame_numbers(:)==frame); % get index in full framelist for current frame
        if ~isempty(index_a)
            x=new_marbles_comingFromRight(marblenum).row_of_centers(index_a); % x left most pixel
            y=new_marbles_comingFromRight(marblenum).col_of_centers(index_a); % y of left most pixel
            matGroundTruthStats(frame,1)=matGroundTruthStats(frame,1)+1;
            %Look for our marbles within a distance of x and y
            index_detected=find([Tracking_marbles.frame]==frame);
            if ~isempty(index_detected)
                x_detected=Tracking_marbles(index_detected).xcenter;
                y_detected=Tracking_marbles(index_detected).ycenter;
                for i_detected=1:length(x_detected)
                    if (iDetectionRadius>=abs(sqrt(((x_detected(i_detected)-x)^2+(y_detected(i_detected)-y)^2))))
                        % If it was detected within 10 pixels:
                        % Count it
                        matGroundTruthStats(frame,2) = matGroundTruthStats(frame,2) + 1;
                        % Add the distance to the sum for averaging
                        sum_of_distances_from_ground_truth = sum_of_distances_from_ground_truth + abs(sqrt(((x_detected(i_detected)-x)^2+(y_detected(i_detected)-y)^2)));
                        break;
                    end
                    
                end
            end
        end
    end
    
    for marblenum=1:size(new_marbles_comingFromLeft,2)       % loop over the individual marbles coming from the left
        index_a=find(new_marbles_comingFromLeft(marblenum).frame_numbers(:)==frame); % get index in full framelist for current frame
        if ~isempty(index_a)
            x=new_marbles_comingFromLeft(marblenum).row_of_centers(index_a); % x left most pixel
            y=new_marbles_comingFromLeft(marblenum).col_of_centers(index_a); % y of left most pixel
            matGroundTruthStats(frame,1)=matGroundTruthStats(frame,1)+1;
            %Look for our marbles within a distance of x and y
            index_detected=find([Tracking_marbles.frame]==frame);
            if ~isempty(index_detected)
                x_detected=Tracking_marbles(index_detected).xcenter;
                y_detected=Tracking_marbles(index_detected).ycenter;
                for i_detected=1:length(x_detected)
                    if (iDetectionRadius>=abs(sqrt(((x_detected(i_detected)-x)^2+(y_detected(i_detected)-y)^2))))
                        % If it was detected within 10 pixels:
                        % Count it
                        matGroundTruthStats(frame,2) = matGroundTruthStats(frame,2) + 1;
                        % Add the distance to the sum for averaging
                        sum_of_distances_from_ground_truth = sum_of_distances_from_ground_truth + abs(sqrt(((x_detected(i_detected)-x)^2+(y_detected(i_detected)-y)^2)));
                        break;
                    end
                    
                end
            end
        end
    end
    
    %% Tracking statistics: count the number of correct pairings in consecutive frames.
    %{
    % For each marble detected in this frame...
    marblesDetectedInThisFrame = marblesWithIDs(frame, :, :);
    for marble = 1 : size(marblesDetectedInThisFrame,2)-1
        if marblesDetectedInThisFrame(1, marble, 1) ~= 0 % Ignore marbles that were not detected in this frame
            id_in_this_frame
    %}
    
end

%% Tracking statistics
for marble = 1 : length(new_marbles_comingFromLeft)
    for appearance = 1 : length(new_marbles_comingFromLeft(marble))
        
        
        
    end
end

for i=1:length(matGroundTruthStats(:,1))
    fprintf('Frame %d. Images in Ground Truth: %d. Detected marbles within range: %d. Ratio %.6f\n',i,matGroundTruthStats(i,1),matGroundTruthStats(i,2),matGroundTruthStats(i,2)/matGroundTruthStats(i,1));
end

fprintf('Final Statistics.\nRatio of detected marbles/ground truth:%.6f\n\n',sum(matGroundTruthStats(:,2))/sum(matGroundTruthStats(:,1)));
success=true;

mean_distances_from_ground_truth = sum_of_distances_from_ground_truth / sum(matGroundTruthStats(i,2));

end

