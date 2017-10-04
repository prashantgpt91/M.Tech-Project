% LearnStructurePara
% learn structure size change statistics
% learning Spos, Mneg, Sneg
% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2008
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function LearnStructurePara()

FeatDir = [GetRootPath 'Gargs\'];
GargParaDir = [GetRootPath 'Training\Parameters\'];

% postive
TrainImgList = [GetRootPath 'Lists\training_pos.txt'];

[TrainFiles1,TrainFiles2] = textread(TrainImgList,'%s%s%*[^\n]');
TrainNum  = length(TrainFiles1);

N = zeros(TrainNum,1);
M = N;
for i=1:TrainNum,
    N(i) = size(load([FeatDir sprintf('%s.yi.txt',TrainFiles1{i})]),1);
    M(i) = size(load([FeatDir sprintf('%s.yi.txt',TrainFiles2{i})]),1);
end

diff = (N-M).*(N-M);
diff = sort(diff);
diff = diff(1:ceil(length(diff)*0.9));   % outlier rejection
Spos = sqrt(sum(diff)/length(diff));

% negative
TrainImgList = [GetRootPath 'Lists\training_neg.txt'];

TrainFiles1 = textread(TrainImgList,'%s%*[^\n]');
TrainNum  = length(TrainFiles1);

N = zeros(TrainNum,1);
for i=1:TrainNum,
    N(i) = size(load([FeatDir sprintf('%s.yi.txt',TrainFiles1{i})]),1);
end

Mneg = mean(N);
Sneg = std(N);

% write
fid = fopen([GargParaDir 'struct.txt'],'w');
fprintf(fid,'%f %f %f',[Spos Mneg Sneg]);
fprintf(fid,'\r\n');
fclose(fid);
