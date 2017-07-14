

function CornerDetection()

[files] = textread([GetRootPath 'lists\list_all.txt'],'%s%*[^\n]');
ImgDir = [GetRootPath 'Images\'];
CornerDir = [GetRootPath 'Corners\'];
margin = 24;    

len = length(files);

for i=1:len,
    disp(['Extracting corners: ' files{i} ' (' num2str(i) '/' num2str(length(files)) ')' ]);
    
    file_name = [ImgDir files{i} '.jpg'];
    
    img = imread(file_name,'jpg');
    
    ss = size(img);
    h = ss(1);  w = ss(2);

    
    if(length(ss)==3) 
        img = rgb2gray(img);
    end
    

    [cim, r, c] = harris(img,3,50,3,0);

    data = [c r];
    %disp('haha')
    %disp(data)

    ind = find(c>margin & c<w-margin & r>margin & r<h-margin);
    
    c = c(ind); 
    r = r(ind);
    data = [r c];
   % disp('hahahuhuhuhu')
   % disp(data)

    fid = fopen([CornerDir files{i} '.txt'],'w');
    for j=1:length(r),    
        fprintf(fid,'%d\t%d\r\n',data(j,1),data(j,2));
    end
    
    fclose(fid);
%    pause;
end
