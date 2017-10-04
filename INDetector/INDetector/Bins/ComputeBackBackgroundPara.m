% ComputeBackBackgroundPara
% Estimate background parameters with negative class training data
% anealing typically is set to 0

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function ComputeBackBackgroundPara(anealing)

TrainImgList = [GetRootPath 'Lists\training_neg.txt'];
FeatDir = [GetRootPath 'Gargs\'];
BackDir = [GetRootPath 'Training\\Parameters\\'];

[FeatDim,N2FeatDim] = GetGlobalParameters;

[TrainFiles] = textread(TrainImgList,'%s%*[^\n]');
TrainNum  = length(TrainFiles);

% 1-node feature files
% read in all features
for i=1:TrainNum,
    FeatFile{i} = load([FeatDir sprintf('%s.yi.txt',TrainFiles{i})]);
end

% compute background parameters 
Y = [];
for i=1:TrainNum,
    Y = [Y ; FeatFile{i}(:,2:(FeatDim+1))];
end

[row,col] = size(Y);
bmean = mean(Y);
%bvar  = var(Y);

% anealed version
bvar = var(Y) + ones(1,FeatDim)*2*anealing/row;

% save background parameters
fid = fopen([BackDir sprintf('V9999.mean.txt')],'w');
fprintf(fid,'%d ',9999);
fprintf(fid,'%f ',bmean);
fprintf(fid,'\r\n');
fclose(fid);
    
fid = fopen([BackDir sprintf('V9999.var.txt')],'w');
fprintf(fid,'%d ',9999);
fprintf(fid,'%f ',bvar);
fprintf(fid,'\r\n');
fclose(fid);

if N2FeatDim~=0
    
    clear Y;

% 2-node feature files
% read in all features
    for i=1:TrainNum,
        FeatFile{i} = load([FeatDir sprintf('%s.yij.txt',TrainFiles{i})]);
    end

% compute background parameters 
    Y = [];
    for i=1:TrainNum,
        Y = [Y ; FeatFile{i}(:,5:(N2FeatDim+4))];
    end
    bmean = mean(Y);
    bvar  = var(Y);

% save background parameters
    fid = fopen([BackDir sprintf('E9999.mean.txt')],'w');
    fprintf(fid,'%d %d ',9999,9999);
    fprintf(fid,'%f ',bmean);
    fprintf(fid,'\r\n');
    fclose(fid);

    fid = fopen([BackDir sprintf('E9999.var.txt')],'w');
    fprintf(fid,'%d %d ',9999,9999);
    fprintf(fid,'%f ',bvar);
    fprintf(fid,'\r\n');
    fclose(fid);

    fid = fopen([BackDir sprintf('E0000.var.txt')],'w');
    fprintf(fid,'%d %d ',9999,9999);
    fprintf(fid,'%f ',bvar);
    fprintf(fid,'\r\n');
    fclose(fid);

end
