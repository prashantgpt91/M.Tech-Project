% GargMatchTest
% Test ARG matching

[dummy,result] = dos('INDetector.exe -m Parameters.txt 19980503_ABC_002578 19980503_CNN_003418');

disp(['match_result = ' result]);
