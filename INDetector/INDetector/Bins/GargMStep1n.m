% RargMStep1n
% Estimate RARG parameters
% RARG learning M-step 1 node
% learning the \mu_i, \sigma_i, and r

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function GargMstep1n(iter)

FeatDir = [GetRootPath 'Gargs\'];
GargParaDir = [GetRootPath 'Training\Parameters\'];
TrainImgList = [GetRootPath 'Lists\training_pos.txt'];
GargXDir = [GetRootPath 'Training\X\'];

[FeatDim,N2FeatDim] = GetGlobalParameters;
[TrainFiles1,TrainFiles2] = textread(TrainImgList,'%s%s%*[^\n]');
TrainNum  = length(TrainFiles1);

% read in all features
for i=1:TrainNum,
    FeatFile1{i} = load([FeatDir sprintf('%s.yi.txt',TrainFiles1{i})]);
    FeatFile2{i} = load([FeatDir sprintf('%s.yi.txt',TrainFiles2{i})]);
end

% read in all one node x
X = [];
x_sum = zeros(TrainNum,1);
% DataWeights = [];
for i=1:TrainNum,
    xs = load([GargXDir sprintf('xiu_%04d.txt',i-1)]);
    [row col] = size(xs);
    x_sum(i) = sum(xs(:,3))/row;
    xs = [ones(row,1)*(i-1) xs];
    X = [X ; xs];
end

% compute parameters
xs  = X;
    
[row col] = size(xs);

Y1 = [];
Y2 = [];

for j=1:row,
   Y1 = [Y1 ; FeatFile1{xs(j,1)+1}(xs(j,2)+1,:)];
   Y2 = [Y2 ; FeatFile2{xs(j,1)+1}(xs(j,3)+1,:)];
end
        
Y1 = Y1(:,2:(FeatDim+1));
Y2 = Y2(:,2:(FeatDim+1));

weight = xs(:,4);
BigW = ones(FeatDim,1)*weight';

Y1 = Y1';
Y2 = Y2';
        
%  mean
	numerator = ((Y2-Y1).*BigW)*(Y2-Y1)';
	denomator = sum(weight);
        
	fvar = numerator/denomator;
	fvar = diag(fvar);
	fid = fopen([GargParaDir 'V0000.var.txt'],'w');
	fprintf(fid,'%d ',0);
	fprintf(fid,'%f ',fvar);
	fprintf(fid,'\r\n');
	fclose(fid);

% Update copy probabilities
aver_x=0;
for i=1:TrainNum,
    aver_x = aver_x+x_sum(i);
end

CopyProb = aver_x/TrainNum;
if iter==1,
    CopyProb = 0.9;
end

fid = fopen([GargParaDir 'r.txt'],'w');
fprintf(fid,'%f ',CopyProb);
fprintf(fid,'\r\n');
fclose(fid);
