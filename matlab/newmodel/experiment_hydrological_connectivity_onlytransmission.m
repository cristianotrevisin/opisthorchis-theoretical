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

% Fish population scaling 0.3
Nf = (A*1000).^0.3;
hf = ones(SC.nNodes,1)*0.001;


%% HYDROLOGICAL CONNECTIVITY
WF = zeros(SC.nNodes,SC.nNodes);
for nn = 1:SC.nNodes
    temp = find(downNode==nn);
    WF(temp,nn) = 1;
    WF(nn,temp) = 1;
end


%% EPIMODEL
% define parameters

par.lambda_F = 0.001;
par.mu_F = 1e-5;

%define setup
setup.nNodes = SC.nNodes;
setup.WF = WF;
setup.H = P';
setup.Nf = Nf;
setup.outlet = zeros(SC.nNodes,1); setup.outlet(A == max(A)) = 1;
setup.Ns = 20000*ones(SC.nNodes,1);


Time = 1:100*365;

%initial solution
y0=zeros(3,SC.nNodes);

y0(3,43) = 100;


Y = model_diffusion(Time,par,setup,y0);

y = Y(:,3:3:end);

%%

figure()
subplot(2,1,1)
imagesc(assignSC(MSC,y(end,:))')
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
plot(Time/365,sum(y,2))

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