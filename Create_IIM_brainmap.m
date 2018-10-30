% WRITTEN BY: ELHUM A SHAMSHIRI (Elhum.Shamshiri@unige.ch)
% CODE SHOULD NOT BE DISTRIBUTED AND ANALYSIS SHOULD NOT BE CONDUCTED WITHOUT PRIOR CONSENT FROM  ELHUM A SHAMSHIRI
% If you would like to use this software for publication please contact

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

clear
clc
close all
%matlabpool;
addpath('D:\Intrisic_Ignition_Model\');
cd E:\FreeSurfer\freesurfer\subjects\subject1\anat\ % INPUT!
%allROIs = dir('*compact.nii');
load('D:\Intrisic_Ignition_Model\names_of_ROIs_reduced_reorg.mat')

for i = 1:length(names_of_ROIs)
    allROIs{i,1} = strcat('roi',num2str(i),'.img');
end

save_fig_lh = 'I:\Data\Patients\subject1\fMRI_Analysis\IIM\lh.fig'; % INPUT!
save_fig_rh = 'I:\Data\Patients\subject1\fMRI_Analysis\IIM\rh.fig'; % INPUT!
load('I:\Data\Patients\subject1\fMRI_Analysis\IIM\IDMI_1sd.mat') % INPUT!
load('I:\Data\Patients\subject1\fMRI_Analysis\IIM\parcellation_matrix_reduced.mat')% INPUT!
%IDMI = single(IDMI); % to make things go faster change from double to single

%% Standardise IDMI values to z-scores
IDMI = standardise_IDMI(IDMI);

%% GET RID OF NAN REGIONS %%
tStart1 = tic;
%orig_num = length(allROIs);
orig_num = size(names_of_ROIs,1);

% row_has_all_zeros = ~any(y, 2);
% a = find(row_has_all_zeros); % find coordinates of zeros in IDMI
% % IDMI(isnan(IDMI)) = []; % get rid of NaNs in IDMI
% for i=1:length(a)
%     allROIs(a(i)-(i-1))=[]; % get rid of NaN regions in allROIs
% end
for i=1:orig_num; % make labels for regions
    labels{i,1} = names_of_ROIs{i,1}(4:end-11);
end
% for i=1:length(a)
%     labels(a(i))=[]; % get rid of NaN regions in labels
% end

%% LOAD DATA %%
for i = 1:length(allROIs)
    imgname = fullfile(pwd,allROIs{i});
    s = load_untouch_nii(imgname);
    img{i} = s.img;
end

%% LEFT HEMISPHERE %%

elliecolour = jet(length(allROIs)); % for the whole brain
figure(1);
[rnk] = floor(tiedrank(IDMI(1:length(IDMI)))); % ranking data order
for i = 1:length(allROIs)*0.5
    brain = double(img{1,i});
    for j = 1:size(brain,1)
        for k = 1:size(brain,2)
            for l = 1:size(brain,3)
                if brain(j,k,l)>0
                    brain(j,k,l)=IDMI(1,i);
                end
            end
        end
    end
    [x,y,z] = ndgrid(1:size(brain,1), 1:size(brain,2), 1:size(brain,3));
    keep = (brain(:) ~= 0); % gets rid of zeros to make outisde of brain transparent
        x = x(keep);
        y = y(keep);
        z = z(keep);
        hplotl(1,i) = plot3(x(:),y(:),z(:),'.','color',elliecolour(rnk(i),:),'MarkerSize',15); % you need to find regions in x,y,z
        hold on
        hl{i} = labels{i,1};
end
hold off
view([-10.5000,68])
colorbar('westoutside');
caxis([min(IDMI) max(IDMI)]);
grid off
box on
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
title('Ignition-Driven Mean Integration (IDMI): Left Hemisphere','fontweight','bold')
gridLegend(hplotl,6,hl,'Location','southoutside');

%saveas(gcf,save_fig_lh)

%% RIGHT HEMISPHERE %%

figure(2);
for i = length(allROIs)*0.5+1:length(allROIs)
    brain = double(img{1,i});
    for j = 1:size(brain,1)
        for k = 1:size(brain,2)
            for l = 1:size(brain,3)
                if brain(j,k,l)>0
                    brain(j,k,l) = IDMI(1,i); % replaces values with IDMI value
                end
            end
        end
    end
    [x,y,z] = ndgrid(1:size(brain,1), 1:size(brain,2), 1:size(brain,3));
    keep = (brain(:) ~= 0); % gets rid of zeros to make outisde of brain transparent
        x = x(keep);
        y = y(keep);
        z = z(keep);
        hplotr(1,i) = plot3(x(:),y(:),z(:),'.','color',elliecolour(rnk(i-(length(allROIs)*0.5)),:),'MarkerSize',15); % you need to find regions in x,y,z
        hold on
        hr{i} = labels{i,1};
end
hold off
view([-170,-74])
colorbar('westoutside');
caxis([min(IDMI) max(IDMI)]);
grid off
box on
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
title('Ignition-Driven Mean Integration (IDMI): Right Hemisphere','fontweight','bold')
hr = hr(length(allROIs)*0.5+1:length(allROIs));
hplotr = hplotr(length(allROIs)*0.5+1:length(allROIs));
gridLegend(hplotr,6,hr,'Location','southoutside');

%saveas(gcf,save_fig_rh)

%% TIMING %%
tEnd1 = toc(tStart1);
fprintf('%d minutes and %f seconds\n',floor(tEnd1/60),rem(tEnd1,60));

%matlabpool close;
