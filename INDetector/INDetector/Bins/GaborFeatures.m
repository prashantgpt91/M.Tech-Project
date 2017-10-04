% Get gabor features at certain location
% The output is a matrix, with each row a feature vector
% corner_id feature_c1 feature_c2 ...

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function vec = GaborFeatures(image,row,col,g_kernel_a,g_kernel_b)

len = length(row);
knum = length(g_kernel_a);
[irow,icol] = size(image);

vec = [];
for i=1:len,
    r = row(i);
    c = col(i);
    
    margin_tag  = 0;
    vecrow = [];
    for k=1:knum,
        w = size(g_kernel_a{k},1);
        rad = (w-1)/2;

        if r-rad<=0 | r+rad>irow | c-rad<=0 | c+rad>icol,
            margin_tag = 1;
            
            aaa
            break;
        end
        
        subimg = image(r-rad:r+rad,c-rad:c+rad);
        a = g_kernel_a{k}(:)'*subimg(:);
        b = g_kernel_b{k}(:)'*subimg(:);
        vecrow = [vecrow sqrt(a*a+b*b)];
    end
    
    if margin_tag>0,
        continue;
    end
    
%    vecrow = [i vecrow];
    vec = [vec; vecrow];
end
