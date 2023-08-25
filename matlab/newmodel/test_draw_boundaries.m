clc
clearvars
close all

load subcatchments

COLS = ["#ff595e", "#ff924c", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"];


CCR = randi(25,[1 max(CTC)]);

CTD = zeros(1,length(CTC));
for i = 1:length(CTD)
    CTD(i) = CCR(CTC(i));
end

CTM = zeros(300,200); MSC = zeros(300,200);
XI = ceil(X);
YI = ceil(Y);

for i = 1:length(X)
    CTM(XI(i),YI(i)) = CTD(i);
    MSC(XI(i),YI(i)) = CTC(i);
end


load FD
load SC
FD.nNodes = length(CTD);
FD.outlet = 19801;
FD.A = FD_A; FD.X = FD_X; FD.Y = FD_Y; FD.downNode = FD_downNode;
% AvailableNodes = setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
thrA = 120; cellsize = 1;



CTA = zeros(300,200);

load AG
for i = 1:length(X)
    CTA(XI(i),YI(i)) = A(CTC(i));
end



figure()
img = image(MSC','CDataMapping','scaled');
colorbar
set(gca,'YDir','normal')
hold on
title('Watersheds')
% for i = 1:FD.nNodes
%     if i ~= FD.outlet
%         if (FD.A(i)>=thrA && ...
%                 abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
%                 abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
%             line([FD.X(i) FD.X(FD.downNode(i))],...
%                 [FD.Y(i) FD.Y(FD.downNode(i))],...
%                 'linewidth', 0.5+4.5*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
%                 'color','white')
%         end
%     end
% end

for ix = 1:size(MSC,1)
    for iy = 2:size(MSC,2)
        if MSC(ix,iy) ~= MSC(ix,iy-1)
            line([ix-0.5 ix+0.5], [iy-0.5 iy-0.5], 'Color','black','linewidth',0.5)
        end
    end
end
for iy = 1:size(MSC,2)
    for ix = 2:size(MSC,1)
        if MSC(ix,iy) ~= MSC(ix-1,iy)
            line([ix-0.5 ix-0.5], [iy-0.5 iy+0.5], 'Color','black','linewidth',0.5)
        end
    end
end