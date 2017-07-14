% GargLearnBatch
% Learning ARG matching 

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

[FeatDim,N2FeatDim] = GetGlobalParameters;

% clear the likelihood file
fid = fopen([GetRootPath 'Results\likelihood.txt'],'w');
fclose(fid);

% set parameter
fid = fopen([GetRootPath 'Training\Parameters\para.txt'],'w');  % Init synchronization file
fprintf(fid,'%d\r\n',FeatDim);
fprintf(fid,'%d\r\n',N2FeatDim);
fclose(fid);

ComputeBackBackgroundPara(0);  % parameters of negative hypothesis 
LearnStructurePara;             % learn size change statistics

disp('Start EM process ...');

% Start E-M process
dos('INDetector.exe -0 Parameters.txt');   % initialization

for i=1:8,   % run 8 iteration
    
% M step
   GargMstep1n(i);  % learn 1n parameters
   
   if(N2FeatDim~=0 && i~=1)
      GargMstep2n(i);  % learn 2n parameters
   end
   
% E step
   dos('INDetector.exe -1 Parameters.txt');
    
   disp(sprintf('Run %d finished\r\n',i));
end
