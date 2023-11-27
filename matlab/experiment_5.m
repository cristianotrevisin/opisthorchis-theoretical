% In here we want to test the effect of fish market by running a simulation


clc; close all; clearvars;
% read OCN

OCN_A = build_OCN("OCN_A.mat",30*10000*10000);
OCN_B = build_OCN("OCN_B.mat",30*10000*10000);
OCN_C = build_OCN("OCN_C.mat",30*10000*10000);

par0 = common_parameters();
par0.lambda_F = 0.01;


[upsA, downsA, labels] = sensitivity_analysis_EE(OCN_A,par0, 3108);
[upsB, downsB, labels_B] = sensitivity_analysis_EE(OCN_B,par0, 2507);
[upsC, downsC, labels_C] = sensitivity_analysis_EE(OCN_C,par0, 0808);

UP = (upsA+upsB+upsC)/3;
DOWN = (downsA + downsB + downsC)/3;

% Get order based on difference
DIFF = abs(UP(1,:)-DOWN(1,:));
[~,IDX] = sort(DIFF,'descend');
UP_TEMP = UP;
DOWN_TEMP = DOWN;
fields_temp = labels;
UP = zeros(4,length(fields_temp),3);
DOWN = zeros(4,length(fields_temp),3);
labels = cell(size(fields_temp));

for rk = 1:length(DIFF)
        UP(:,rk,1) = upsA(:,IDX(rk));
        UP(:,rk,2) = upsB(:,IDX(rk));
        UP(:,rk,3) = upsC(:,IDX(rk));

        DOWN(:,rk,1) = downsA(:,IDX(rk));
        DOWN(:,rk,2) = downsB(:,IDX(rk));
        DOWN(:,rk,3) = downsC(:,IDX(rk));
        labels(rk) = fields_temp(IDX(rk));
end


%%
lower = min(min(UP,[],'all'),min(DOWN,[],'all'))*100;
upper = max(max(UP,[],'all'),max(DOWN,[],'all'))*100;
absolute = max(abs(lower),abs(upper));
figure
tiledlayout(2,2)
for i = 1:4
    nexttile
    d = barh(1:length(labels),squeeze(DOWN(i,:,:))*100);
    hold on
    for j = 1:3; d(j).FaceColor = 'blue'; d(j).EdgeColor = 'none'; end;
    u = barh(1:length(labels),squeeze(UP(i,:,:))*100);
    for j = 1:3; u(j).FaceColor = 'red'; u(j).EdgeColor = 'none'; end;
    xlabel('[%] variation')
    box off
    ax = gca;
    ax.YColor = 'w';
    ax.YAxis.Label.Color='k';
    yticks(1:length(labels))
    yticklabels(labels)
    set(gca,'FontSize',8)
    set(gca,'YDir','reverse')
    if i == 1
        title('Worm burden in humans')
    elseif i == 2
        title ('Egg concentration')
    elseif i == 3
        title ('Prevalence of infected snails')
    else
        title('Cyst burden in fish')
    end
    xlim([-absolute, absolute])
end
