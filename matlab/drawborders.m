function drawborders(MSC,BorderColor)
if nargin == 1
    BorderColor = 'black';
end
hold on
for ix = 1:size(MSC,1)
    for iy = 2:size(MSC,2)
        if MSC(ix,iy) ~= MSC(ix,iy-1)
            line([ix-0.5 ix+0.5], [iy-0.5 iy-0.5], 'Color',BorderColor,'linewidth',0.5)
        end
    end
end
for iy = 1:size(MSC,2)
    for ix = 2:size(MSC,1)
        if MSC(ix,iy) ~= MSC(ix-1,iy)
            line([ix-0.5 ix-0.5], [iy-0.5 iy+0.5], 'Color',BorderColor,'linewidth',0.5)
        end
    end
end

end