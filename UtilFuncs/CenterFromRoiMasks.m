function []=CenterFromRoiMasks(roiMask,cell_ids,opt,arg1, arg2)
%Renamed from CenterFromGlomLabel.m
%Input is roiMask / Counting mask whose value is 0-n and non-zero
%components form patches
%This function calculate the center of mass of each patch and numbering in
%in figures
if ~exist('arg1','var')
arg1='r';
end
if ~exist('arg2','var')
arg2=8;
end   
text_c= arg1;
text_size = arg2;


num=max(roiMask(:));

for i=1:num
    clear pos
   [y,x]=find(roiMask==i);
   
   pos_x(i)=mean(x);
   pos_y(i)=mean(y);
end

if opt~=0
    for i=1:num
        text(pos_x(i),pos_y(i),num2str(cell_ids(i)),'FontSize',text_size,'Color',text_c,'FontWeight','bold')
    end
end