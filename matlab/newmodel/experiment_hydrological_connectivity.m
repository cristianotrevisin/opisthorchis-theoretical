clc
clearvars
close all


load ../../dataOCN/subcatchments

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


rng(2021)



load ../../dataOCN/FD
load ../../dataOCN/SC
FD.nNodes = length(CTD);
FD.outlet = 1921;
FD.A = FD_A; FD.X = FD_X; FD.Y = FD_Y; FD.downNode = FD_downNode;
% AvailableNodes = setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
thrA = 120; cellsize = 1;



CTA = zeros(NX,NY);

load ../../dataOCN/AG
for i = 1:length(X)
    CTA(XI(i),YI(i)) = A(CTC(i));
end


% DEFINING POPULATION
% population to be defined by sorting areas wrt accumulation flow
% [~,ranks]=ismember(A,sort(A,'descend'));

% population is sorted randomly
ranks = randperm(length(A));
% Populations are Zipf distributed
P = zipf(ranks, 2, 250);

%to downstream acc
[~,IDSORT] = sort(A);
P(IDSORT) = sort(P);


%%
M = zeros(SC.nNodes,SC.nNodes);
DST = zeros(SC.nNodes,SC.nNodes);
LOC = 1*ones(1,SC.nNodes);
par.phi = 1;

for sc1 = 1:SC.nNodes
    for sc2 = 1:SC.nNodes
        DST(sc1,sc2) = sqrt((SCX(sc1)-SCX(sc2))^2+(SCY(sc1)-SCY(sc2))^2);
        if sc1~=sc2
            M(sc1,sc2) = P(sc1)*P(sc2)/(DST(sc1,sc2)^par.phi);
        end
    end
end
M = M./(sum(M,1)).*(1-LOC);
for sc = 1:SC.nNodes; M(sc,sc) = LOC(sc); end
% figure()
% imagesc(M)



%% HYDROLOGICAL CONNECTIVITY
%WS = zeros(SC.nNodes,SC.nNodes);
W = zeros(SC.nNodes,SC.nNodes);
for nn = 1:SC.nNodes
    temp = find(downNode==nn);
    W(temp,nn) = 1;
    %WF(temp,nn) = 1;
    %WF(nn,temp) = 1;
end
figure(); imagesc(W)

%% EPIMODEL
% define parameters
par.beta_FH = 4.898e-2; 
par.beta_HS = 9.160e-2;
par.beta_SF = 3.477e-5*(A*1000).^0.3./mean((A*1000).^0.3);
par.alpha = 1;

par.mu_H = 1/365/10;
par.mu_S = 1/365;
par.mu_F = 1/365/2.5;

par.lambda_FU = 0.0001;
par.lambda_FD = 0.0002;
par.mu_F = 1e-5;

%define setup
setup.nNodes = SC.nNodes;
setup.Cf = 0.001*ones(SC.nNodes,1);
setup.outlet = zeros(SC.nNodes,1); setup.outlet(A == max(A)) = 1;

setup.W = W;
setup.M = M;
setup.H = P';
setup.N_F = 1000*ones(SC.nNodes,1);

Time = 1:50*365;

%initial solution
y0=zeros(3,SC.nNodes);

y0(1:3,43) = [0; 0; 1];



y = model(Time,par,setup,y0);

WH = y(:,1:3:end);
ES = y(:,2:3:end);
CF = y(:,3:3:end);

%%
figure()
subplot(2,1,1)
imagesc(assignSC(MSC,CF(end,:))')
set(gca,'YDir','normal')
drawborders(MSC)
colorbar
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA && ...
                abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
                abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
            line([FD.X(i) FD.X(FD.downNode(i))],...
                [FD.Y(i) FD.Y(FD.downNode(i))],...
                'linewidth', 0.5+4.5*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
                'color','white')
        end
    end
end


subplot(2,1,2)
plot(Time/365,sum(CF.*setup.N_F',2))

%%

figure()
subplot(3,1,1)
imagesc(assignSC(MSC,WH(end,:))')
set(gca,'YDir','normal')
drawborders(MSC)
colorbar
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA && ...
                abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
                abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
            line([FD.X(i) FD.X(FD.downNode(i))],...
                [FD.Y(i) FD.Y(FD.downNode(i))],...
                'linewidth', 0.5+4.5*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
                'color','white')
        end
    end
end
title('Worms in humans')

subplot(3,1,2)
imagesc(assignSC(MSC,ES(end,:))')
set(gca,'YDir','normal')
drawborders(MSC)
colorbar
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA && ...
                abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
                abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
            line([FD.X(i) FD.X(FD.downNode(i))],...
                [FD.Y(i) FD.Y(FD.downNode(i))],...
                'linewidth', 0.5+4.5*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
                'color','white')
        end
    end
end
title('Eggs in snails')

subplot(3,1,3)
imagesc(assignSC(MSC,CF(end,:))')
set(gca,'YDir','normal')
drawborders(MSC)
colorbar
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA && ...
                abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
                abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
            line([FD.X(i) FD.X(FD.downNode(i))],...
                [FD.Y(i) FD.Y(FD.downNode(i))],...
                'linewidth', 0.5+4.5*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
                'color','white')
        end
    end
end
title('Cercariae in fish')

CFS = CF; CFS(CF<0) = -1; CFS(CF>0) = 1; 
figure()
imagesc(assignSC(MSC,CFS(end,:))')
set(gca,'YDir','normal')
drawborders(MSC)
colorbar
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA && ...
                abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
                abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
            line([FD.X(i) FD.X(FD.downNode(i))],...
                [FD.Y(i) FD.Y(FD.downNode(i))],...
                'linewidth', 0.5+4.5*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
                'color','white')
        end
    end
end
title('Cercariae in fish')

%%


% shading interp
% figure()
% surf(xx,yy,WHP','edgecolor','none')
%  view(0,90)
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