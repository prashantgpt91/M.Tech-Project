% Estimate RARG parameters
% RARG learning M-step 2 node

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function RargMstep2n(iter)

sampling = 1000;

FeatDir = [GetRootPath 'Gargs\\'];
GargParaDir = [GetRootPath 'Training\\Parameters\\'];
TrainImgList = [GetRootPath 'Lists\\training_pos.txt'];
GargXDir = [GetRootPath 'Training\\X\\'];

[FeatDim,N2FeatDim] = GetGlobalParameters;
[TrainFiles1,TrainFiles2] = textread(TrainImgList,'%s%s%*[^\n]');
TrainNum  = length(TrainFiles1);

% read in all features
for i=1:TrainNum,
    FeatFile1{i} = load([FeatDir sprintf('%s.yij.txt',TrainFiles1{i})]);
    FeatFile2{i} = load([FeatDir sprintf('%s.yij.txt',TrainFiles2{i})]);
end

% read in all two node x
X = [];
for i=1:TrainNum,
    xs = load([GargXDir sprintf('xijuv_%04d.txt',i-1)]);
    [row col] = size(xs);
    xs = [ones(row,1)*(i-1) xs];
    X = [X ; xs];
end

xs = X;
% compute parameters
non_neg = find(xs(:,2)~=-1);
xs = xs(non_neg,:);
non_neg = find(xs(:,4)~=-1);
xs = xs(non_neg,:);

[row col] = size(xs);

sample_id = randperm(row);
sample_id = sample_id(1:sampling);
xs = xs(sample_id,:);
row = sampling;

Y1 = [];
Y2 = [];
for j=1:row,
    
    dir = xs(j,3);
    if dir<0,
        disp('error in X, dir1 <0');
    end
    
    dir = xs(j,5);
    if(dir>0)
        n1 = xs(j,8);
        n2 = xs(j,9);
    else
        n1 = xs(j,9);
        n2 = xs(j,8);
    end
        
    if n1>=n2,
        disp('error in X output');
    end
        
    edge_id1 = xs(j,2);
    edge_id2 = xs(j,4);
    yf1 = FeatFile1{xs(j,1)+1};
    yf2 = FeatFile2{xs(j,1)+1};
        
    Y1 = [Y1 ; yf1(edge_id1+1,5:6)*dir]; 
    Y2 = [Y2 ; yf2(edge_id2+1,5:6)*dir];
end

weight = xs(:,10);
BigW = ones(N2FeatDim,1)*weight';

Y1 = Y1';
Y2 = Y2';
    
%   mean
numerator = ((Y2-Y1).*BigW)*(Y2-Y1)';
denomator = sum(weight);

%   numerator = Ysub*diag(weight)*Ysub';
fvar  = numerator/denomator;
fvar  = diag(fvar);

fid = fopen([GargParaDir 'E0000.var.txt'],'w');
fprintf(fid,'%d %d ',0,0);
fprintf(fid,'%f ',fvar);
fprintf(fid,'\r\n');
fclose(fid);
