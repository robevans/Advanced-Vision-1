
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
[matMarbles,success]=detect_marbles(directory);

if (~success)
    return;
end

[Tracking_marbles,success]=track_marbles(matMarbles);

if (~success)
   return;
end   

success=display_tracking(Tracking_marbles,directory);

if (~success)
    return;
end    

success=ground_truth_statistics(directory,Tracking_marbles)


end

