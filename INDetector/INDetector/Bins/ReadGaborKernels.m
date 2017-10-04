% read gabor kernels from dir

% Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
% Contact : 
%			Dong-Qing Zhang :	 dongqing@gmail.com, 
%			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu

function [g_kernel_a,g_kernel_b] = ReadGaborKernels()

KernelDir = '.\GaborKernels\';
scal_no = 3;
orien_no = 4;

n = 1;

for i=1:scal_no,
    for j=1:orien_no,
        file_name = [KernelDir 'gabor_' sprintf('%02d',i-1) '_' sprintf('%02d',j-1) '_a.txt'];
        g_kernel_a{n} = load(file_name);

        file_name = [KernelDir 'gabor_' sprintf('%02d',i-1) '_' sprintf('%02d',j-1) '_b.txt'];
        g_kernel_b{n} = load(file_name);
        
        n = n+1;
    end
end
