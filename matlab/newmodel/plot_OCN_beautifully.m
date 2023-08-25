clc
clearvars
close all

rng(2021)

load ../../dataOCN/subcatchments

load ../../dataOCN/outlet

COLS = ["#ff595e", "#ff924c", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"];

NX = ceil(max(X));
NY = ceil(max(Y));

CCR = randi(25,[1 max(CTC)]);

CTD = zeros(1,length(CTC));
for i = 1:length(CTD)
    CTD(i) = CCR(CTC(i));
end

CTM = zeros(NX,NY); MSC = zeros(NX,NY); ELE = zeros(NX,NY);
XI = ceil(X);
YI = ceil(Y);

load ../../dataOCN/FD




% fare bei disegni


load ../../dataOCN/SC
FD.nNodes = length(CTD);
FD.outlet = outlet;
FD.A = FD_A; FD.X = FD_X; FD.Y = FD_Y; FD.downNode = FD_downNode; FD.Z = FD_Z;



% Snail population
% Suitability 

b1 = 0.8;
b2 = 0.005;

Rice_Pond = zeros(size(FD.Z));
for i = 1:length(FD.Z)
    if rand <= b1*exp(-b2*FD.Z(i))
        Rice_Pond(i) = 1;
    end
end



CTA = zeros(NX,NY); Rice_Pond_MAT= zeros(NX,NY);

load ../../dataOCN/AG
for i = 1:length(X)
    CTA(XI(i),YI(i)) = A(CTC(i));
end

% Fish population scaling 0.3
Nf = (A*1000).^0.3;
hf = ones(SC.nNodes,1)*0.001;

thrA = 55; cellsize = 1;
for i = 1:length(X)
    CTM(XI(i),YI(i)) = CTD(i);
    MSC(XI(i),YI(i)) = CTC(i);
    ELE(XI(i),YI(i)) = FD.Z(i);
    Rice_Pond_MAT(XI(i),YI(i)) = Rice_Pond(i);
end





h=figure();
%img = imagesc(-ELE','CDataMapping','direct','AlphaData',0.4);
img = imagesc(CTM','CDataMapping','direct','AlphaData',0.4);

drawborders(MSC,'k')
set(gca,'YDir','normal')
hold on

colriver = hex2rgb('#2276ea');
colriver(4) = 0.5;

lineWidth = 0.6*(FD.A*0.1).^0.4;
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA)
            line([FD.X(i) FD.X(FD.downNode(i))],[FD.Y(i) FD.Y(FD.downNode(i))],...
                    'color','#2776ea','linewidth',0.9*(FD.A(i)*1e-1)^0.3,'AlignVertexCenters','on')
        elseif (FD.A(i)>=2) && (FD.A(i)<thrA)
            line([FD.X(i) FD.X(FD.downNode(i))],[FD.Y(i) FD.Y(FD.downNode(i))],...
                    'color',colriver,'linewidth',0.5*(FD.A(i)*1e-1)^0.3,'AlignVertexCenters','on')
        end
    end
end
line([0.5 0.5],[0.5 size(MSC,2)+.5],'color','k')
line([size(MSC,1)+.5 size(MSC,1)+.5],[0.5 size(MSC,2)+.5],'color','k')
line([0.5 size(MSC,1)+.5],[0.5 0.5],'color','k')
line([0.5 size(MSC,1)+.5],[size(MSC,2)+.5 size(MSC,2)+.5],'color','k')
set(gcf,'GraphicsSmoothing','on')
axis equal
axis off

% DRAW POPULATIONS
% plot 2 - human population randomly distributed
% plot 3 - downstream accumulation for human population
% plot 4 - fish concentration; increasing downstream (plot relationship)
% plot 1 - catchment number

% population is sorted randomly
ranks = randperm(length(A));
% Populations are Zipf distributed
P = zipf(ranks, 2, 250);

%to downstream acc
[~,IDSORT] = sort(A);
PDA = zeros(size(P));
PDA(IDSORT) = sort(P);








figure()
colormap lines
img = imagesc(CTM','CDataMapping','scaled','AlphaData',0.4);
drawborders(MSC,'k')
set(gca,'YDir','normal')
hold on
for sc = 1:SC.nNodes
    text(SCX(sc),SCY(sc),num2str(sc),'Color','black','FontName','Fedra Sans Pro','FontSize',9)
end
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA)
            line([FD.X(i) FD.X(FD.downNode(i))],[FD.Y(i) FD.Y(FD.downNode(i))],...
                    'color','#2776ea','linewidth',0.9*(FD.A(i)*1e-1)^0.3,'AlignVertexCenters','on')
        end
    end
end
line([0.5 0.5],[0.5 size(MSC,2)+.5],'color','k')
line([size(MSC,1)+.5 size(MSC,1)+.5],[0.5 size(MSC,2)+.5],'color','k')
line([0.5 size(MSC,1)+.5],[0.5 0.5],'color','k')
line([0.5 size(MSC,1)+.5],[size(MSC,2)+.5 size(MSC,2)+.5],'color','k')
set(gcf,'GraphicsSmoothing','on')
axis equal
axis off

%%
[XXX,YYY] = meshgrid(unique(FD.X),unique(FD.Y));
figure()

surf(XXX',YYY',ELE/50)
shading interp
hold on
patch([XXX(1,end); XXX(1,1); XXX(1,:)'],[YYY(1,end); YYY(1,1); YYY(1,:)'],[0; 0; ELE(:,1)/50],'k')
patch([XXX(end,1); XXX(1,1); XXX(:,1)],[YYY(end,1); YYY(1,1); YYY(:,1)],[0 0 ELE(1,:)/50],'k')
axis equal
axis off

%%
%%
figure()
title('Snail suitability (random)')
colormap sky
img = imagesc(Rice_Pond_MAT','CDataMapping','direct','AlphaData',0.4);
drawborders(MSC,'white')
set(gca,'YDir','normal')
hold on
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA)
            line([FD.X(i) FD.X(FD.downNode(i))],[FD.Y(i) FD.Y(FD.downNode(i))],...
                    'color','#2776ea','linewidth',0.9*(FD.A(i)*1e-1)^0.3,'AlignVertexCenters','on')
        end
    end
end
line([0.5 0.5],[0.5 size(MSC,2)+.5],'color','k')
line([size(MSC,1)+.5 size(MSC,1)+.5],[0.5 size(MSC,2)+.5],'color','k')
line([0.5 size(MSC,1)+.5],[0.5 0.5],'color','k')
line([0.5 size(MSC,1)+.5],[size(MSC,2)+.5 size(MSC,2)+.5],'color','k')
axis equal
axis off

%%

figure()
title('Human population (random)')
colormap sky
imagesc(assignSC(MSC,log10(P))','CDataMapping','scaled','AlphaData',0.4);
colorbar
drawborders(MSC,'white')
set(gca,'YDir','normal')
hold on
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA)
            line([FD.X(i) FD.X(FD.downNode(i))],[FD.Y(i) FD.Y(FD.downNode(i))],...
                    'color','#2776ea','linewidth',0.9*(FD.A(i)*1e-1)^0.3,'AlignVertexCenters','on')
        end
    end
end
line([0.5 0.5],[0.5 size(MSC,2)+.5],'color','k')
line([size(MSC,1)+.5 size(MSC,1)+.5],[0.5 size(MSC,2)+.5],'color','k')
line([0.5 size(MSC,1)+.5],[0.5 0.5],'color','k')
line([0.5 size(MSC,1)+.5],[size(MSC,2)+.5 size(MSC,2)+.5],'color','k')
axis equal
axis off

figure()
title('Human population (downstream accumulation)')
colormap sky
imagesc(assignSC(MSC,log10(PDA))','CDataMapping','scaled','AlphaData',0.4);
colorbar
drawborders(MSC,'white')
set(gca,'YDir','normal')
hold on
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA)
            line([FD.X(i) FD.X(FD.downNode(i))],[FD.Y(i) FD.Y(FD.downNode(i))],...
                    'color','#2776ea','linewidth',0.9*(FD.A(i)*1e-1)^0.3,'AlignVertexCenters','on')
        end
    end
end
line([0.5 0.5],[0.5 size(MSC,2)+.5],'color','k')
line([size(MSC,1)+.5 size(MSC,1)+.5],[0.5 size(MSC,2)+.5],'color','k')
line([0.5 size(MSC,1)+.5],[0.5 0.5],'color','k')
line([0.5 size(MSC,1)+.5],[size(MSC,2)+.5 size(MSC,2)+.5],'color','k')
axis equal
axis off

figure()
title('Fish concentration')
colormap parula
imagesc(assignSC(MSC,Nf)','CDataMapping','scaled','AlphaData',0.4);
set(gca,'YDir','normal')
drawborders(MSC,'white')
hold on
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA)
            line([FD.X(i) FD.X(FD.downNode(i))],[FD.Y(i) FD.Y(FD.downNode(i))],...
                    'color','#2776ea','linewidth',0.9*(FD.A(i)*1e-1)^0.3,'AlignVertexCenters','on')
        end
    end
end
line([0.5 0.5],[0.5 size(MSC,2)+.5],'color','k')
line([size(MSC,1)+.5 size(MSC,1)+.5],[0.5 size(MSC,2)+.5],'color','k')
line([0.5 size(MSC,1)+.5],[0.5 0.5],'color','k')
line([0.5 size(MSC,1)+.5],[size(MSC,2)+.5 size(MSC,2)+.5],'color','k')
axis equal
axis off


%%
function M = assignSC(MSC,X)
    M = zeros(size(MSC));
    for iz = 1:length(X)
        M(MSC==iz) = X(iz);
    end
end


function P = zipf(rank, expn, minP)
%     ranks = 1:1:N;
%     pmf = (ranks.^(-expn))/sum(ranks.^(-expn));
%     samples = rand(1,M);
%     p = cumsum(pmf(:));
%     [~,x] = histc(samples,[0;p/p(end)])

    H = sum(1./(rank.^expn));
    
    
    P = 1./rank.^expn./H;
    P = P/min(P)*minP;
    
end
