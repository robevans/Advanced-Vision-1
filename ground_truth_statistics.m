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
            x_detected=Tracking_marbles{frame}(:,1);
            y_detected=Tracking_marbles{frame}(:,2);
            %Look for our marbles within a distance of x and y
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
%            end
        end
    end
    
    for marblenum=1:size(new_marbles_comingFromLeft,2)       % loop over the individual marbles coming from the left
        index_a=find(new_marbles_comingFromLeft(marblenum).frame_numbers(:)==frame); % get index in full framelist for current frame
        if ~isempty(index_a)
            x=new_marbles_comingFromLeft(marblenum).row_of_centers(index_a); % x left most pixel
            y=new_marbles_comingFromLeft(marblenum).col_of_centers(index_a); % y of left most pixel
            matGroundTruthStats(frame,1)=matGroundTruthStats(frame,1)+1;
            x_detected=Tracking_marbles{frame}(:,1);
            y_detected=Tracking_marbles{frame}(:,2);
            %Look for our marbles within a distance of x and y
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
%            end
        end
    end
    
end

n_correct_pairings = 0;
n_erroneous_pairings = 0;
first_assigned_label = 0;
second_assigned_label = 0;

%% Tracking statistics: Count the number of correct and erroneous pairings out of the marbles that we detected.
for marbleSet = 1 : 2
    if marbleSet == 1
        ground_truth_marbles = new_marbles_comingFromLeft;
    else
        ground_truth_marbles = new_marbles_comingFromRight;
    end
    for marble = 1 : length(ground_truth_marbles)
            for appearance = 1 : length(ground_truth_marbles(marble).frame_numbers)-1
                % Ground truth information for the marble in CURRENT frame
                gt_x = ground_truth_marbles(marble).row_of_centers(appearance);
                gt_y = ground_truth_marbles(marble).col_of_centers(appearance);
                frame = ground_truth_marbles(marble).frame_numbers(appearance);
                
                % Check if we detected that marble within 10 pixels in CURRENT frame
                detectedMarbles = marblesWithIDs(frame,:,:);
                for ith_detected = 1 : size(detectedMarbles, 2)
                    x_detected = detectedMarbles(1, ith_detected, 1);
                    y_detected = detectedMarbles(1, ith_detected, 2);
                    if ~(x_detected == 0 && y_detected == 0)
                        if (iDetectionRadius >= abs(sqrt(((x_detected-gt_x)^2 + (y_detected-gt_y)^2))))
                            % If we did detect it within 10 pixels:
                            first_assigned_label = ith_detected;
                        end
                    end
                end
                
                % Ground truth information for the same marble in the NEXT frame
                gt_x = ground_truth_marbles(marble).row_of_centers(appearance+1);
                gt_y = ground_truth_marbles(marble).col_of_centers(appearance+1);
                frame = ground_truth_marbles(marble).frame_numbers(appearance+1);
                
                % Check if we detected that marble within 10 pixels in the NEXT frame
                detectedMarbles = marblesWithIDs(frame,:,:);
                for ith_detected = 1 : size(detectedMarbles, 2)
                    x_detected = detectedMarbles(1, ith_detected, 1);
                    y_detected = detectedMarbles(1, ith_detected, 2);
                    if ~(x_detected == 0 && y_detected == 0)
                        if (iDetectionRadius >= abs(sqrt(((x_detected-gt_x)^2 + (y_detected-gt_y)^2))))
                            % If we did detect it within 10 pixels:
                            second_assigned_label = ith_detected;
                        end
                    end
                end
                
                % If we gave the marble identical labels in consecutive
                % frames then it was tracked properly.
                if first_assigned_label == second_assigned_label
                    % We have a correct pairing by the tracker!
                    n_correct_pairings = n_correct_pairings + 1;
                else
                    % Oh no!
                    n_erroneous_pairings = n_erroneous_pairings + 1;
                end
                
            end
    end
end

fprintf('\nDetection Statistics:\n  Number of marbles in ground truth: %i\n  Ratio of detected marbles/ground truth: %i/%i\n\n',sum(matGroundTruthStats(:,1)),sum(matGroundTruthStats(:,2)),sum(matGroundTruthStats(:,1)));

mean_distances_from_ground_truth = sum_of_distances_from_ground_truth / sum(matGroundTruthStats(:,2))

fprintf('Tracking Statistics:\n  Percentage of detected marbles that were correctly tracked in consecutive frames: %.2f%%\n  Number of correct pairings: %i\n  Number of erroneous pairings: %i\n\n',n_correct_pairings/(n_correct_pairings+n_erroneous_pairings)*100,n_correct_pairings,n_erroneous_pairings);
success=true;
end

