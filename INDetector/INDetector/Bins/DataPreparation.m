 % DataPreparation

% CreateDirectories;

CornerDetection;   % Extract interest points using harris corner detector

FeatureExtraction;  % Extract features around interest points 

FeaturesToGargs_txt(6,0);  % bag-of-parts representation, 6 is the minum size of an ARG or BoP
% FeaturesToGargs_txt(6,1);  % ARG representation
