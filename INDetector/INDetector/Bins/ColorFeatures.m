% Get averaged color features at the specified location
% The output is a matrix, with each row a feature vector
% corner_id feature_c1 feature_c2 ...

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function vec = ColorFeatures(image,row,col)

rad = 24;  % eqauls to cut-out margin for corner detection, see CornerDetection.m
len = length(row);
[irow,icol] = size(image);

vec = [];
for i=1:len,
    r = row(i);
    c = col(i);
    
    vecrow = [];

    if r-rad<=0 | r+rad>irow | c-rad<=0 | c+rad>icol,
        Disp('The cut-out margin has to be larger than 24, see CornerDetection.m');
        continue;
    end

%     imshow(image), hold on;
%     plot(c,r,'r+');
%     hold off;
    
    subimg = image(r-rad:r+rad,c-rad:c+rad,:);
    r = subimg(:,:,1);
    g = subimg(:,:,2);
    b = subimg(:,:,3);
    
    vecrow = [mean(r(:)) mean(g(:)) mean(b(:))];
    
%    vecrow = [i vecrow];
    vec = [vec; vecrow];
end
