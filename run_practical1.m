
function [ success ] = run_practical1( directory )
%AV Practical 1
%January 30, 2014       Creation
%function run_practical1
%   This function runs the Advanced Vision practical 1
%Input : directory - the directory name where images are
%Output : success - success state for last step executed
%Authors
%Robert Evans
%Francisco Aguirre

%[Images_marbles,success]=detect_marbles_adaptive(directory);
[matMarbles,vecnDetected,success]=detect_marbles(directory,false);

if (~success)
    return;
end
%TODO remove after completion
local_matMarbles=matMarbles;
local_directory=directory;
local_vecnDetected=vecnDetected;
[Tracking_marbles,matMarblesPosition,nTrackedMarbles,success]=track_marbles(matMarbles,directory, vecnDetected);

local_Tracking_marbles=Tracking_marbles;
local_matMarbles_position=matMarblesPosition;
if (~success)
   return;
end   

success=display_tracking(directory,nTrackedMarbles,matMarblesPosition);

if (~success)
    return;
end    

success=ground_truth_statistics(directory,Tracking_marbles);
if (~success)
    fprintf('Ground Truth statistics could not be calculated.');
    return;
end    


end

