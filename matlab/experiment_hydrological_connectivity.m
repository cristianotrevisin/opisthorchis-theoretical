clc
clearvars
close all


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


rng(2021)

% ergodicity
% modello spazialmente implicito 
% iverse repliche del modello SE 
% effetto della pesca
% vedere effeto della distribuzione random delaal popolazione contro
% concentrazione downstreamn
% ipotesi What if 
% espericmenti numerici con diverse strutture della popolazione
% effetto della pressione sulla pesca fishing effort 
% modello di trasport e-g. tempo di percrrenza del pesce percato fino al
% mercaot 
% effetto di diverse comdizioni e.g. temperatura piogge, etc 
% confronto tra esplicito e implicito 


% fare bei disegni

load ../dataOCN/FD
load ../dataOCN/SC
FD.nNodes = length(CTD);
FD.outlet = 1921;
FD.A = FD_A; FD.X = FD_X; FD.Y = FD_Y; FD.downNode = FD_downNode;
% AvailableNodes = setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
thrA = 120; cellsize = 1;



CTA = zeros(NX,NY);

load ../dataOCN/AG
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

% Fish population scaling 0.3
Nf = (A*1000).^0.3;
hf = ones(SC.nNodes,1)*0.001;

% figure()
% img = image(CTM','CDataMapping','scaled');
% set(gca,'YDir','normal')
% hold on
% title('Watersheds')
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
% 
% for sc = 1:SC.nNodes
%     plot(SCX(sc),SCY(sc),'.r','MarkerSize',0.5+1.5*log(P(sc)))
%     text(SCX(sc),SCY(sc),num2str(sc),'Color','white')
% end
% 
% set(gca,'visible','off')
% xlim([0 NX])
% ylim([0 NY])

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
WS = zeros(SC.nNodes,SC.nNodes);
WF = zeros(SC.nNodes,SC.nNodes);
for nn = 1:SC.nNodes
    temp = find(downNode==nn);
    WS(nn,temp) = 1;
    WF(temp,nn) = 1;
    WF(nn,temp) = 1;
end


%% EPIMODEL
% define parameters
par.beta_FH = 0;%4.898e-5; 
par.beta_HS = 0;%9.160e-11;
par.beta_SF = 0;%3.477e-5;

par.mu_H = 0;%1/365/10;
par.mu_S = 0;%1/365;
par.mu_F = 0;%1/365/2.5;

par.lambda_S = 0;
par.lambda_F = 0.01;

%define setup
setup.nNodes = SC.nNodes;
setup.Cf = hf;
setup.WS = WS;
setup.WF = WF;
setup.M = M;
setup.H = P';
setup.Nf = Nf;
setup.Ns = 20000*ones(SC.nNodes,1);



Time = 1:400*365;

%initial solution
y0=zeros(3,SC.nNodes);

y0(1:3,43) = [0; 0; .5];
%y0(1,253)=0; %i
y = model(Time,par,setup,y0);

WH = y(:,1:3:end);
SI = y(:,2:3:end);
FI = y(:,3:3:end);

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
title('Worm burden in humans')

subplot(3,1,2)
imagesc(assignSC(MSC,SI(end,:))')
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
title('Quota of infected snails')

subplot(3,1,3)
imagesc(assignSC(MSC,FI(end,:))')
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
title('Quota of infected fish')

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