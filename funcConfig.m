function [ ParameterValue ] = funcConfig( sParameterName )
%function funcConfig
%This function receives a parameter name and returns its value
%Useful to keep configuration variables in one place

switch sParameterName
    case 'nImagesAverageBackground'
        ParameterValue = 10;
    case 'sDirectory'
        ParameterValue = './SEQ1/';
    case 'fileAdaptiveBackground'
        ParameterValue = 'fileAdaptiveBackground.mat';
    case 'nMaxMarbles'
        ParameterValue = 18;
    case 'nNumHypothesis'
        ParameterValue=100;
end;

end
