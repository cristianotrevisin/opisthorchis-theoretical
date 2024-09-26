% SENSITIVITY ANALYSIS

clc; close all; clearvars;
% read OCN
OCN_A = build_OCN("OCN_A.mat",30*10000*10000);
OCN_B = build_OCN("OCN_B.mat",30*10000*10000);
OCN_C = build_OCN("OCN_C.mat",30*10000*10000);

par0 = common_parameters();


[upsA, downsA, labels] = sensitivity_analysis_EE(OCN_A,par0, 3108);
[upsB, downsB, labels_B] = sensitivity_analysis_EE(OCN_B,par0, 2507);
[upsC, downsC, labels_C] = sensitivity_analysis_EE(OCN_C,par0, 0808);

UP = (upsA+upsB+upsC)/3;
DOWN = (downsA + downsB + downsC)/3;

% Get order based on difference
DIFF = abs(UP(1,:,1,1)-DOWN(1,:,1,1));
[~,IDX] = sort(DIFF,'descend');
UP_TEMP = UP;
DOWN_TEMP = DOWN;
fields_temp = labels;
UP = zeros(3,length(fields_temp),3,3,2);
DOWN = zeros(3,length(fields_temp),3,3,2);
labels = cell(size(fields_temp));

for rk = 1:length(DIFF)
        UP(:,rk,1,:,:) = upsA(:,IDX(rk),:,:);
        UP(:,rk,2,:,:) = upsB(:,IDX(rk),:,:);
        UP(:,rk,3,:,:) = upsC(:,IDX(rk),:,:);

        DOWN(:,rk,1,:,:) = downsA(:,IDX(rk),:,:);
        DOWN(:,rk,2,:,:) = downsB(:,IDX(rk),:,:);
        DOWN(:,rk,3,:,:) = downsC(:,IDX(rk),:,:);
        labels(rk) = fields_temp(IDX(rk));
end

save("Sensitivity_Analysis_RES.mat")

%% GENERATE FIGURES
UP(:,end-2:end,:,:,:) = [];
DOWN(:,end-2:end,:,:,:) = [];
labels(end-2:end) = [];
lower = min(min(UP,[],'all'),min(DOWN,[],'all'))*100;
upper = max(max(UP,[],'all'),max(DOWN,[],'all'))*100;
absolute = max(abs(lower),abs(upper));
figure
tiledlayout(3,2)
for i = 1:3
    for m = 1:2
        nexttile
        d20 = barh(1:length(labels),squeeze(DOWN(i,:,:,1,m))*100);
        hold on
        d10 = barh(1:length(labels),squeeze(DOWN(i,:,:,2,m))*100);
        d05 = barh(1:length(labels),squeeze(DOWN(i,:,:,3,m))*100);
        for j = 1:3; d20(j).FaceColor = '#033270'; d20(j).EdgeColor = 'none'; end;
        for j = 1:3; d10(j).FaceColor = '#1368aa'; d10(j).EdgeColor = 'none'; end;
        for j = 1:3; d05(j).FaceColor = '#4091c9'; d05(j).EdgeColor = 'none'; end;
        u20 = barh(1:length(labels),squeeze(UP(i,:,:,1,m))*100);
        u10 = barh(1:length(labels),squeeze(UP(i,:,:,2,m))*100);
        u05 = barh(1:length(labels),squeeze(UP(i,:,:,3,m))*100);
        for j = 1:3; u20(j).FaceColor = '#65010c'; u20(j).EdgeColor = 'none'; end;
        for j = 1:3; u10(j).FaceColor = '#cb1b16'; u10(j).EdgeColor = 'none'; end;
        for j = 1:3; u05(j).FaceColor = '#f26a4f'; u05(j).EdgeColor = 'none'; end;
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
            title ('Prevalence of infected snails')
        else
            title('Cyst burden in fish')
        end
        xlim([-absolute, absolute])
        set(gca,'FontSize',9)
    end
end
