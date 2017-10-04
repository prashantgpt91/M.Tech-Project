% Extract Features using corners
% spatial+color+gabor, 2+3+12=17 dimension

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function FeatureExtraction()

[files] = textread([GetRootPath 'lists\list_all.txt'],'%s%*[^\n]');
ImgDir = [GetRootPath 'Images\'];
CornerDir = [GetRootPath 'Corners\'];
FeatureDir = [GetRootPath 'Features\'];

[g_kernel_a,g_kernel_b] = ReadGaborKernels;

len = length(files);

for i=1:len,

        disp(['Extracting features: ' files{i} ' (' num2str(i) '/' num2str(length(files)) ')' ]);
    
% read image
        file_name = [ImgDir files{i} '.jpg'];

        img = imread(file_name,'jpg');

        ss = size(img);
        h = ss(1);  w = ss(2);
         
 % read corners
        file_name = [CornerDir files{i} '.txt'];
        Corners = load(file_name);

        point_list_size = size(Corners,1);
        
        if point_list_size~=0,
           point_list = Corners(:,1:2);
           
% spatial features
           svec = [point_list(:,2) h-point_list(:,1)];
           
% color features            
           cvec = ColorFeatures(img,point_list(:,1),point_list(:,2));
    
           if(length(ss)==3)  % color image
                img = rgb2gray(img);
           end
           img = double(img)/255;
           
% gabor features        
           point_list = Corners(:,1:2);
           gvec = GaborFeatures(img,point_list(:,1),point_list(:,2),g_kernel_a,g_kernel_b);
            
        end
        
       fid = fopen([FeatureDir files{i} '.txt'],'w');
       
       if point_list_size~=0,
            fprintf(fid,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\r\n',[(1:point_list_size)' svec cvec gvec]');
       end  
       
       fclose(fid);

end
