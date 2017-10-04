% CornerDetection
% this script is used for generating corners

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function CornerDetection()

[files] = textread([GetRootPath 'lists\list_all.txt'],'%s%*[^\n]');
ImgDir = [GetRootPath 'Images\'];
CornerDir = [GetRootPath 'Corners\'];
margin = 24;     % margin for garbor filter kernel, see line 41 and ColorFeatures.m

len = length(files);

for i=1:len,
    disp(['Extracting corners: ' files{i} ' (' num2str(i) '/' num2str(length(files)) ')' ]);
    
    file_name = [ImgDir files{i} '.jpg'];
    
    img = imread(file_name,'jpg');
    
    ss = size(img);
    h = ss(1);  w = ss(2);

%    imshow(img), hold on;
    
    if(length(ss)==3)  % color image
        img = rgb2gray(img);
    end
    
% please download harris.m from http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/index.html    
    [cim, r, c] = harris(img,4,500,4,0);

    data = [c r];
    
%     if files{i}(10)~='C',
%         ind = find(c>24 & c<w-24 & r>24 & r<h-24);  % for TREC ABC news
%     else
%         ind = find(c>24 & c<w-24 & r>24 & r<h-45);  % remove TREC CNN news ticker
%     end
    ind = find(c>margin & c<w-margin & r>margin & r<h-margin);
    
    c = c(ind); 
    r = r(ind);
    data = [r c];

%    plot(c,r,'r+'), title('corners detected');
%        hold off;
    
    fid = fopen([CornerDir files{i} '.txt'],'w');
    for j=1:length(r),    
        fprintf(fid,'%d\t%d\r\n',data(j,1),data(j,2));
    end
    
    fclose(fid);
%    pause;
end
