clc
clearvars
close all

rng(2021)

load ../dataOCN/subcatchments

COLS = ["#ff595e", "#ff924c", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"];

NX = ceil(max(X));
NY = ceil(max(Y));

CCR = randi(25,[1 max(CTC)]);

CTD = zeros(1,length(CTC));
for i = 1:length(CTD)
    CTD(i) = CCR(CTC(i));
end

CTM = zeros(NX,NY); MSC = zeros(NX,NY);
XI = ceil(X);
YI = ceil(Y);

for i = 1:length(X)
    CTM(XI(i),YI(i)) = CTD(i);
    MSC(XI(i),YI(i)) = CTC(i);
end


% fare bei disegni

load ../dataOCN/FD
load ../dataOCN/SC
FD.nNodes = length(CTD);
FD.outlet = 1921;
FD.A = FD_A; FD.X = FD_X; FD.Y = FD_Y; FD.downNode = FD_downNode;

thrA = 10; cellsize = 1;



CTA = zeros(NX,NY);

load ../dataOCN/AG
for i = 1:length(X)
    CTA(XI(i),YI(i)) = A(CTC(i));
end


h=figure();
img = image(CTM','CDataMapping','direct','AlphaData',0.2);
colormap lines
set(gca,'YDir','normal')
hold on
title('Watersheds')
lineWidth = 0.9*(FD.A*0.1).^0.4;
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA)

            line([FD.X(i) FD.X(FD.downNode(i))],[FD.Y(i) FD.Y(FD.downNode(i))],...
                    'color','b','linewidth',0.9*(FD.A(i)*1e-1)^0.4,'AlignVertexCenters','on')
        end
    end
end
set(gcf,'GraphicsSmoothing','on')
