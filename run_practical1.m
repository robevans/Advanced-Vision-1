%This function runs the practical 1
function [ success ] = run_practical1( directory )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%AV Practical 1
%Robert Evans
%Francisco Aguirre

[Images_marbles,success]=detect_marbles(directory);

if (~success)
    return;
end

[Tracking_marbles,success]=track_marbles(Images_marbles);

if (~success)
   return;
end   

success=display_tracking(Tracking_marbles,directory);

if (~success)
    return;
end    

success=ground_truth_statistics(Tracking_marbles)


end

