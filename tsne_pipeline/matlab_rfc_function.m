
function [feature_bar,feature_scatter] = matlab_rfc_function(in_table,palette)
%in table is a table of the correct format
%palette is a desired color format

in_cluster = cellfun(@(x) num2str(x),num2cell(in_table.Cluster),'UniformOutput', false);
in_data = in_table{:,2:end};

predictor_names = in_table.Properties.VariableNames(2:end);
b = TreeBagger(100,in_data,in_cluster,'PredictorNames',predictor_names,'OOBVarImp','On','OOBPredictorImportance','on','NumPredictorsToSample','all','PredictorSelection','interaction-curvature');
% figure
% plot(oobError(b))
%
% xlabel('Number of Grown Trees')
% ylabel('Out-of-Bag Classification Error')

font = 'Arial';
feature_bar= figure(33);
set(feature_bar,'units','normalized','outerposition',[0.1 0.1 .5 .5],'DefaultTextFontName', font, 'DefaultAxesFontName',font)
[ordered_imp,I] = sort(b.OOBPermutedPredictorDeltaError,'descend');
sorted_names = predictor_names(I);

bar(ordered_imp)

ylabel('Feature Importance')
h = gca;
h.XTick = 1:length(predictor_names);
h.XTickLabel = sorted_names;
h.XTickLabelRotation = 90;
%h.FontSize = 24;
h.FontWeight = 'bold';

h.XColor = 'k';
h.YColor = 'k';
h.Box = 'off';
h.LineWidth = 2;
ytick_old = h.YTick;
%h.YLim(2) = h.YTick(end)+(h.YTick(end)-h.YTick(end-1))*.25;
h.YTick = ytick_old;
h.TickDir = 'out';
title('Idenficication of Important Features')
% 
feature_scatter = figure(22);
hold on
[s,e] =mdsProx(fillProximities(b),'Colors','rbgk');
clusters = unique(in_table.Cluster);
disp(clusters)
for k =1:length(clusters')
    sub_table = in_table(in_table.Cluster == clusters(k),:);
    scatter(sub_table.(sorted_names{1}),sub_table.(sorted_names{2}),'filled','MarkerEdgeColor',palette(k,:),'MarkerFaceColor',palette(k,:))
    hold on
end

xlabel(['Normed, Centered ',sorted_names{1}], 'FontWeight','bold','Color','k')
ylabel(['Normed, Centered ',sorted_names{2}], 'FontWeight','bold','Color','k')
legend(cellfun(@(x,y)[x,y], repmat({'Cluster '},1,k), cellfun(@(x) {num2str(x)},num2cell(1:k)) ,'UniformOutput',false))

 ax_scatter= findall(feature_scatter,'Type','axes')

%set(ax_scatter,'fontsize',12)
%set(ax_scatter,'fontweight','bold')
set(ax_scatter,'GridColor','k')

ax_scatter.XLim(1) = ax_scatter.XLim(1)*1.1;
ax_scatter.YLim(1) = ax_scatter.YLim(1)*1.1;
ax_scatter.XLim(2) = ax_scatter.XLim(2)*1.1;
ax_scatter.YLim(2) = ax_scatter.YLim(2)*1.1;
ax_scatter.XColor = 'k';
ax_scatter.YColor = 'k';
ax_scatter.LineWidth = 2;
ax_scatter.TickDir = 'out';


title('Scatter: Important Features');

%fig_yer.PaperPosition = [1 1 10 10 ];
%%print(fullfile('.','results_fao_goetzman','RFC_MATLAB_OUTPUT',all_diseases{i}),'-depsc');

end

