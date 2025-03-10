function [intense_mask, intense_mask_nand] = IntensityMask(profiles,threshold,op)

if nargin==2
    op = 10;
end
[Nx, Ny, Npc] = size(profiles);
complex_sum = abs(sum(profiles,3));

complex_sum_feature = StandardizeImage(complex_sum);

intense_mask = 1*(complex_sum_feature>threshold);
intense_mask = bwareaopen(intense_mask,op);


intense_mask_nand = 1*intense_mask;

intense_mask_nand(intense_mask==0) = nan;

end




function data = StandardizeImage(data)
sizerarr = size(data);
data = data(:);
data = data-mean(data);
data = data/(std(data));
data = reshape(data, sizerarr);
end