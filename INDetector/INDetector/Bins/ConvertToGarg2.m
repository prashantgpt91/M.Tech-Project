% convert features to gargs
% need to change if use different features
% called by FeaturesToGargs_txt

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function garg_out = ConvertToGarg2(features)

[row,col] = size(features);
if row==0,
    garg_out.n = [];
    garg_out.e = [];
    return;
end

garg_out.n = [features(:,1) features(:,2:16)];
%garg_out.n = features;

edge_arr = zeros(row*(row-1)/2,1+1+2+2);

eid = 1;
for i=1:row,
    for j=(i+1):row,
        edge_arr(eid,:) = [eid-1 1 i-1 j-1 features(i,2:3)-features(j,2:3)];
        eid = eid+1;
    end
end
garg_out.e = edge_arr;
