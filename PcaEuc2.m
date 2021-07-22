function [EcuDist, EcuDistAll, EcuDistVcords, EcudDistVcordsMean]= PcaEuc2(Scores, ExpVar, ColorGrp, MarkerGrp, Comp, NetCon, BetweenGrps, RepEcu, Textex, Save)
%
% Versatile representation of PCA scores
%
% Function can utilize data coming form PCA calculations and creates
% 2D PCA with mean Euclidean Distance between groups and within groups.
%
% e.g.:
% PcaEuc2(PCA.SCORE, PCA.EXPLAINED, Type, Batch, [1,2],'yes','yes','yes','Impact of oil origin','yes');
%
%
% OUTPUTS:
%   1. Graphical output ->  Gives PCA with spider/web for group data points
%   2. Numerical output -> Matrix with mean Euclidean distance between
%                           data points coming from that same group
%
%   EcuDist            - Mean Euclidean Distance within group
%   EcuDistAll         - Euclidean Distance within group to mean of group
%   EcuDistVcords      - Euclidean Distance between mean coordinates and 
%                        mean of groups
%   EcudDistVcordsMean - Mean Euclidean Distance coordinates between groups
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
% Scores      - PCA scores data matrix from Matlab's PCA function. 
%               (It could be also different source of scores, 
%               however, they should be prepared in order as in matlab function)
%
% ExpVar      - Column vector or Matlab's PCA function with numerical values
%               of explained variance. Sorted in descending order.
%
% ColorGrp    - Information about sample main grouping. It will be
%               represented as different colors.
%
% MarkerGrp   - Information about sample subgrouping. 
%               It will be represented as different marker. 
%
% Comp        - Matrix with number of components as X,Y axis e.g. [1,2], 
%               which indicates use of PC1 and PC2 on X and Y axis.
%
% NetCon      - 'Yes' or 1 for representation of Euclidean distances
%
% BetweenGrps - 'Yes' or 1 if between groups represent and calculate mean
%               Euclidean distance
%
% RepEcu      - 'Yes' or 1 for representation of Euclidean distances mean
%               in annotation box
%
% Textex      - title of plot and name to save plot
%
% Save        - Input 'Yes' if you want save graphical output into *.TIF file
%
% Author: Wojciech Wojtowicz
% Created: 15.01.2021
% Modified: 01.07.2021
% Contact: wojciech.wojtowicz@pwr.edu.pl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cheaking if input is table
if istable(ColorGrp)
    ColorGrp = table2array(ColorGrp); 
    MarkerGrp = table2array(MarkerGrp);
end

% Checking unique labels in first and secondary groups information
[GrpUniq.One,GrpUniq.Two,GrpUniq.Three] = unique(ColorGrp);
[T1,~,T2] = unique(GrpUniq.Three); 

% Preapring coloring scheme according to number of main groups
cmapLen = length(T1);
cmap = hsv2rgb([(0:cmapLen-1)'/cmapLen,ones(cmapLen,1),.75*ones(cmapLen,1)]);
cmap = flip(cmap);

% Calcaulation of mean cooardinates for each unique group from 1st
for i=1:size(T1)
V(i,:) = [mean(Scores(GrpUniq.Three==T1(i),Comp(1))),mean(Scores(GrpUniq.Three==T1(i),Comp(2)))];
end

% Check RepEcu input style
if isempty(NetCon)
        NetCon = 0;
elseif isstring(NetCon) || iscell(NetCon) || ischar(NetCon)
    if strcmp(NetCon, 'yes') || strcmp(NetCon, 'Yes')
        NetCon = 1;
    else
        NetCon = 0;
    end
elseif isnumeric(NetCon)
    if NetCon ==1
        NetCon = 1;
    else
        NetCon = 0;
    end
end

% Ploting lines between mean coordinates and observations
figure('units','normalized','outerposition',[0 0 1 1]);
hold on
    for j = 1:size(T1)
        EcuDist(j) = mean(sqrt((V(j,1)-Scores(GrpUniq.Three==T1(j),Comp(1))).^2 + (V(j,2)-Scores(GrpUniq.Three==T1(j),Comp(2))).^2));    
        EcuDistAll(:,j) = sqrt((V(j,1)-Scores(GrpUniq.Three==T1(j),Comp(1))).^2 + (V(j,2)-Scores(GrpUniq.Three==T1(j),Comp(2))).^2);  
        if NetCon == 1
            for i = 1:size(Scores(:,Comp(1))) 
            if T2(i) == T1(j)
                cmapInt(j,:) = cmap(j,:)+cmap(j,:)/5;
                plot([V(j,1); Scores(i,Comp(1))], [V(j,2); Scores(i,Comp(2))],'--', 'Color', cmapInt(j,:), 'LineWidth', 0.05,'HandleVisibility','off');
            end
            end
        end
    end

% Ploting mean data points
for i=1:size(T1)
    if NetCon == 1
scatter(mean(Scores(GrpUniq.Three==T1(i),Comp(1))), mean(Scores(GrpUniq.Three==T1(i),Comp(2))), 46, cmap(i,:), 'filled','d','HandleVisibility','off');
    end
V(i,:) = [mean(Scores(GrpUniq.Three==T1(i),Comp(1))),mean(Scores(GrpUniq.Three==T1(i),Comp(2)))];
end

% Mean coordinates for all observations
Vcoord = [mean(Scores(:,Comp(1))), mean(Scores(:,Comp(2)))];

% Ecu. dist. between mean coords and mean of groups from 2nd information
for j=1:size(T1)
EcuDistVcords(j) = mean(sqrt((Vcoord(1,1)-V(j,1)).^2 + (Vcoord(1,2)-V(j,2)).^2));
end

% Mean between ecu. dist.  
EcudDistVcordsMean = mean(EcuDistVcords);

% Check BetweenGrps input style
if isempty(BetweenGrps)
        BetweenGrps = 0;
elseif isstring(BetweenGrps) || iscell(BetweenGrps) || ischar(BetweenGrps)
    if strcmp(BetweenGrps, 'yes') || strcmp(BetweenGrps, 'Yes')
        BetweenGrps = 1;
    else
        BetweenGrps = 0;
    end
elseif isnumeric(BetweenGrps)
    if BetweenGrps ==1
        BetweenGrps = 1;
    else
        BetweenGrps = 0;
    end
end

% Ploting between groups data points and lines between them
if BetweenGrps==1 && NetCon == 1
for j=1:size(T1)
plot(Vcoord(1,1),Vcoord(1,2),'*k','MarkerFaceColor','k');
plot([Vcoord(1,1),V(j,1)], [Vcoord(1,2), V(j,2)],'--k');
end
end

% Placing of annotation box
dim = [0.15 .505 .4 .4];
MnStr= [];

% Strings for annotation box
% Check RepEcu input style
if isempty(RepEcu)
        RepEcu = 0;
elseif isstring(RepEcu) || iscell(RepEcu) || ischar(RepEcu)
    if strcmp(RepEcu, 'yes') || strcmp(RepEcu, 'Yes')
        RepEcu = 1;
    else
        RepEcu = 0;
    end
elseif isnumeric(RepEcu)
    if RepEcu ==1
        RepEcu = 1;
    else
        RepEcu = 0;
    end
end

if RepEcu ==1
for i =1:size(GrpUniq.One)
 if ~iscell(GrpUniq.One(i))
str = ['Mean(Eucl. Dist.) ', num2str(GrpUniq.One(i)) '= ' num2str(round(EcuDist(i),2))];
MnStr = [MnStr; {str}];
 else
str = ['Mean(Eucl. Dist.) ', GrpUniq.One(i) '= ' num2str(round(EcuDist(i),2))];
MnStr = [MnStr; {[str{:}]}];
 end
end
if BetweenGrps==1
str2 = ['Mean(Eucl. Dist.) between= ', num2str(round(EcudDistVcordsMean,2))];
MnStr = [MnStr; str2];
end
annotation('textbox',dim,'String',MnStr,'FitBoxToText','on','FontSize', 20);
end

% Ploting data points according to color from 1st information and marker
% type according to 2nd information
if ~isempty(MarkerGrp)
[GrpUniq1.One,GrpUniq1.Two,GrpUniq1.Three] = unique(MarkerGrp);
[T1_2,~,T2_2] = unique(GrpUniq1.Three); 
if length(T1_2)>length(T1)
for j =1:size(T1_2)
for i=1:size(Scores(:,1),1)
    if T2(i)==T1_2(j)
        GrpCmap(i,:) = cmap(j,:);
    end 
end
end
else
for j =1:size(T1)
for i=1:size(Scores(:,1),1)      
    if T2(i)==T1(j)    
         GrpCmap(i,:) = cmap(j,:);
    end   
end
end
end  
% Finite number of markers - maximum groups in 2nd information
IndMarker =  {'o','s','d','x','v','*','^','+','>','<'};
for j =1:size(T1_2)
for i=1:size(Scores(:,1),1)
     if T2_2(i)==T1_2(j)
         hold on
         if length(IndMarker) >= j
            plot(Scores(i,Comp(1)), Scores(i,Comp(2)),IndMarker{j},'Color', GrpCmap(i,:),'MarkerFaceColor',GrpCmap(i,:), 'MarkerSize',16,'HandleVisibility','off');
         else
             OvsJ=j-10;
             plot(Scores(i,Comp(1)), Scores(i,Comp(2)),IndMarker{OvsJ},'Color', GrpCmap(i,:),'MarkerFaceColor',GrpCmap(i,:), 'MarkerSize',16,'HandleVisibility','off');
         end
     end
end
end
else
[~, ind2] = sort(ColorGrp);
gscatter(Scores(ind2,Comp(1)), Scores(ind2,Comp(2)), ColorGrp(ind2), cmap, [], 45);
end

% Preparing custom legend
if ~isempty(MarkerGrp)
for i = 1:size(T1)
    if iscell(GrpUniq.One)
    Pts = string(GrpUniq.One(i));
    else
     Pts = num2str(GrpUniq.One(i));
    end
h(i) = plot(NaN,NaN,'^','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:),'MarkerSize',16, 'DisplayName', Pts);
clear Pts
end

k =1;
for i = size(T1,1)+1:(size(T1_2)+size(T1,1))
    if iscell(GrpUniq1.One)
        Pts1 = string(GrpUniq1.One(k));
    else
        Pts1 = num2str(GrpUniq1.One(k));
    end
    if length(IndMarker) >= k
        h(i) = plot(NaN,NaN,IndMarker{k},'Color','k','MarkerFaceColor',[1 1 1],'MarkerSize',16, 'DisplayName', Pts1);
        clear Pts1
        k=k+1;
    else
        k=k-10;
        h(i) = plot(NaN,NaN,IndMarker{k},'Color','k','MarkerFaceColor',[1 1 1],'MarkerSize',16, 'DisplayName', Pts1);
        clear Pts1
        k=k+1;       
    end
end
% Figure parameters and saving options
title(Textex,'FontWeight','normal');
xlabel(['PC',num2str(Comp(1)),' ', num2str(round(ExpVar(1,1),2)),' %'], 'FontSize', 20);
ylabel(['PC',num2str(Comp(2)),' ', num2str(round(ExpVar(2,1),2)),' %'], 'FontSize', 20);
box on, grid on;
set(gca, 'FontSize', 24);
alpha 0.5;

plot([0,0], ylim, '--k','HandleVisibility','off');
plot(xlim, [0,0], '--k','HandleVisibility','off');

if ~isempty(MarkerGrp)
legend(h,'Location', 'northeast');
end
switch Save
    case 'yes'
    print([ 'PCA-' Textex],'-dtiff','-r0');
    case 'Yes'
    print([ 'PCA-' Textex],'-dtiff','-r0');
end
end
