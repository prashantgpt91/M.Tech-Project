% GargPadding
% padding the ARG so that its size is at least vertex_num
% called by FeaturesToGargs_txt

function garg_out = GargPaddingMem(garg,vertex_num)

[row,col] = size(garg);

if row==0,
    garg_out = [];
    return;
end
    
if row>=vertex_num,
     garg_out = garg;
     return;
end
    
n = row+1;
for j=1:(vertex_num-row),
    indarr = randperm(row);
    ind = indarr(1);
    garg = [garg; [n garg(ind,2:col)]];
    n = n+1;
end

garg_out = garg;
