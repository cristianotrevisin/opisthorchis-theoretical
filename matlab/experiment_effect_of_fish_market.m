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
thrA = 30; cellsize = 1;



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

figure()
img = image(CTM','CDataMapping','scaled');
set(gca,'YDir','normal')
hold on
title('Watersheds')
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

for sc = 1:SC.nNodes
    plot(SCX(sc),SCY(sc),'.r','MarkerSize',0.5+1.5*log(P(sc)))
    text(SCX(sc),SCY(sc),num2str(sc),'Color','white')
end

set(gca,'visible','off')
xlim([0 NX])
ylim([0 NY])

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
par.beta_FH = 4.898e-5; 
par.beta_HS = 9.160e-11;
par.beta_SF = 3.477e-5;

par.mu_H = 1/365/10;
par.mu_S = 1/365;
par.mu_F = 1/365/2.5;

par.lambda_S = 0;
par.lambda_F = 0.;

%define setup
setup.nNodes = SC.nNodes;
setup.Cf = hf;
setup.WS = WS;
setup.WF = WF;
setup.M = M;
setup.H = P';
setup.Nf = Nf;
setup.Ns = 20000*ones(SC.nNodes,1);

%define equivalent implicit network
setup_I.H = sum(P);
setup_I.Nf = sum(Nf);
setup_I.Ns = sum(setup.Ns);
setup_I.M = 1;
setup_I.Cf = 1;
setup_I.WS = 0;
setup_I.WF = 0;
setup_I.nNodes = 1;

Time = 1:100*365;

%initial solution
y0=zeros(3,SC.nNodes);
%[a, b, c] = EE_OPI(par,setup,3);
y0(1:3,43) = [1/P(43); 0; 0];
%y0(1,253)=0; %initial infected in node 3

%initial solution for implicit
y0_I=zeros(3,1);
%[a, b, c] = EE_OPI(par,setup,3);
y0_I(1:3,1) = [1/sum(P); 0; 0];
%y0(1,253)=0; %initial infected in node 3

y = model(Time,par,setup,y0);
y_I = model(Time,par,setup_I,y0_I);

WH = y(:,1:3:end);
SI = y(:,2:3:end);
FI = y(:,3:3:end);
%compute_indices(par,setup,y,Time);
%%
figure()
plot(Time/365,y_I(:,1))
hold on
plot(Time/365, WH*P'/(sum(P)))
legend('SI','SE')
xlabel('Year')
ylabel('W_H')
%%
lim=100*365;
figure()
yyaxis left
plot(Time(1:lim)/365,y_I(1:lim,1))
ylabel('W_H')
yyaxis right
hold on
plot(Time(1:lim)/365,y_I(1:lim,2))
plot(Time(1:lim)/365,y_I(1:lim,3))
ylabel('i_S or i_F')
xlabel('Year')
legend('WH','IS','IF')
[WI,SII,FII] = EE_OPI(par,setup_I,1)

%%

% find endemic equilibrium for all nodes

EE_vec = zeros(SC.nNodes,1);

for i = 1:SC.nNodes
    I = EE_OPI(par,setup,i);
    if isnan(I) == false
        EE_vec(i) = 1;
    end
end

EE_pixel = assignSC(MSC,EE_vec);

figure()
map = [0.9569    0.9451    0.8706
    0.8784    0.4784    0.3725];
colormap(map)
image(EE_pixel','CDataMapping','scaled')
set(gca,'YDir','normal')
hold on
title('Endemic equilibrium')
drawborders(MSC,'black')
for sc = 1:SC.nNodes
    %plot(SCX(sc),SCY(sc),'.r','MarkerSize',0.5+1.5*log(P(sc)))
    text(SCX(sc),SCY(sc),num2str(sc),'Color','black','FontName','Fedra Sans Pro','FontSize',16)
end
lineWidth = 0.5 + 2 * sqrt(FD.A / (FD.nNodes * cellsize^2));
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA && ...
                abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
                abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
            % line([FD.X(i) FD.X(FD.downNode(i))],...
            %     [FD.Y(i) FD.Y(FD.downNode(i))],...
            %     'linewidth', 0.5+8*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
            %     'color','#0073cf')
             plot([FD.X(i) FD.X(FD.downNode(i))],...
         [FD.Y(i) FD.Y(FD.downNode(i))],...
         'LineWidth', lineWidth(i),...
         'Color', '#0073cf');
        end
    end
end
%%
% population
Nthr = par.mu_F*par.mu_H*par.mu_S/(par.beta_FH*par.beta_HS*par.beta_SF)./setup.Nf./setup.Ns;

figure()
imagesc(assignSC(MSC,Nthr)')
set(gca,'YDir','normal')
drawborders(MSC)
colorbar
title('Threshold population')
lineWidth = 0.5 + 2 * sqrt(FD.A / (FD.nNodes * cellsize^2));
for i = 1:FD.nNodes
    if i ~= FD.outlet
        if (FD.A(i)>=thrA && ...
                abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
                abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
            % line([FD.X(i) FD.X(FD.downNode(i))],...
            %     [FD.Y(i) FD.Y(FD.downNode(i))],...
            %     'linewidth', 0.5+8*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
            %     'color','#0073cf')
             plot([FD.X(i) FD.X(FD.downNode(i))],...
         [FD.Y(i) FD.Y(FD.downNode(i))],...
         'LineWidth', lineWidth(i),...
         'Color', '#0073cf');
        end
    end
end


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
to_plot = [A WH'];
to_plot = sortrows(to_plot,1);

WHP = to_plot(:,2:end);
WHP(WHP==0)=NaN;
WHP = fillmissing(WHP,'linear');

[xx,yy] = meshgrid(sort(A),Time);

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