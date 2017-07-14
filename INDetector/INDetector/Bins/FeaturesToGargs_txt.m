% Convert the single feature file to the text Garg file
% Input : VertexNum -- the minimum vertex number of ARG

function FeaturesToGargs_txt(VertexNum,gen_edge_feature)

drive = GetRootPath;
[files] = textread([drive 'Lists\list_all.txt'],'%s%*[^\n]');
Num = length(files);

n = 0;
for i=1:Num,
    
     disp(['Generating ARGs or BoPs: ' files{i} ' (' num2str(i) '/' num2str(length(files)) ')' ]);
    
     FullPath = [drive 'features\' files{i} '.txt'];
     if ~exist(FullPath),
         files{i}
         n=n+1;
         
         'no such file'
         continue;
     end
     
     garg = load(FullPath);
     garg_out = GargPaddingMem(garg,VertexNum);
     garg_out = ConvertToGarg2(garg_out); 
     
     file1n = [drive 'Gargs\' files{i} '.yi.txt'];
     fid =fopen(file1n,'w');
     len = size(garg_out.n,1);
 % you may need to change the dimension here
     fprintf(fid,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\r\n',garg_out.n');
     fclose(fid);

     if(gen_edge_feature==1)
         file1n = [drive 'Gargs\' files{i} '.yij.txt'];
         fid =fopen(file1n,'w');
         len = size(garg_out.e,1);
         fprintf(fid,'%d %d %d %d %f %f\r\n',garg_out.e');
         fclose(fid);
     end
     
end
