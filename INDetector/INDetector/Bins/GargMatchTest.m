% GargMatchTest
%lambda=0.894
% Test ARG matching

[dummy,result] = dos('INDetector.exe -m Parameters.txt 19980220_ABC_003073 19980430_CNN_040778');
disp(['match_result = ' result]);
