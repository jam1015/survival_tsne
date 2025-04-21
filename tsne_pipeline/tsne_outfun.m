function stop = tsne_outfun(optimValues,state,disease_table,num_clust,cluster_criterion)
% disp(optimValues)
% disp(state)
stop = false; % do not stop by default
color_matrix = [ 31 118 175  ; 255 186 126 ]./255
groups = strrep(disease_table.sample_type_samples,'_',' ');

switch state
    case 'init'
        h = figure;
    case 'iter'
        
        gscatter(optimValues.Y(:,1),optimValues.Y(:,2),groups,color_matrix,'..',[10 10],'off')
        title('Embedding');
       
        drawnow
        
    case 'done'
        
        
        color_map_string = 'aeb01'; %setting color map
        
        palette = GetMatSurvColorPalette(color_map_string);
        
        X = optimValues.Y;  %getting the tsne values
       
      
      
        
        eva = evalclusters(X,'kmeans',cluster_criterion,'KList',num_clust); %seeing how many clusters there are
        [cluster,centroid] = kmeans(optimValues.Y,eva.OptimalK,'EmptyAction','singleton');  %identifying the clusters
        close all
        
        
        
        %% Plotting the clusters
        fh_tsne = figure(1);
        num_clusters = max(max(cluster));
        for k = 1:num_clusters
            scatter(X(cluster==k,1),X(cluster==k,2),'MarkerEdgeColor',palette(k,:),'MarkerFaceColor',palette(k,:))
            
            hold on
        end
        xlabel('tSNE 1')
        ylabel('tSNE 2')
        
        title([strrep(disease_table.disease_type{1},'_',' '),': tSNE Embedding'])
        num_strings = cellfun(@(z) {num2str(z)},num2cell(1:num_clusters));
        cluster_strings = repmat({'Cluster '},1,num_clusters);
        legend(cellfun( @(x,y) {[x,y]}, cluster_strings,num_strings))
        
        
        
        
        %% Performing Survival Analysis on the clusters
        [p,fh_logrank,stats]=MatSurv_tsne(disease_table.x_TIME_TO_EVENT,disease_table.x_EVENT,cellstr([repmat('Cluster ',size(cluster)), num2str(cluster)]),'Title','Cluster Survival','NoRiskTable',true, ...
            'PairwiseP',true,'Print',true,'LineColor',color_map_string,'TimeMin',1,'TimeUnit','Days');
        
        
        %% Using random forest to identify features of the clusters
        
         used_vars = disease_table.Properties.VariableNames;
        gen_var_start = find( strcmp('samples', used_vars))+1;
        gen_var_end   = find( strcmp('x_EVENT', used_vars))-1;
        
        rfc_table  = [table(cluster,'VariableNames',{'Cluster'}),disease_table(:,gen_var_start:gen_var_end)];
        [fh_feature_bar, fh_feature_scatter]=matlab_rfc_function(rfc_table,palette);
       
        
        %% Collecting output plots into one figure
        ax_tsne= findall(fh_tsne,'Type','axes');
        ax_logrank = findall(fh_logrank,'Type','axes');
        ax_feature_bar = findall(fh_feature_bar,'Type','axes');
        ax_feature_scatter = findall(fh_feature_scatter,'Type', 'axes');
        
        fh_output = figure;
        sp1 = subplot(2,2,1);
        sp2 = subplot(2,2,2);
        sp3 = subplot(2,2,3);
        sp4 = subplot(2,2,4);
        % semicolons here change the code behavior. probably has something
        % to do with the output. 

        ax_tsne.Position = sp1.Position
       ax_tsne.OuterPosition = sp1.OuterPosition
        
        ax_logrank.Position      = sp2.Position
        ax_logrank.OuterPosition = sp2.OuterPosition
        
        ax_feature_bar.Position      = sp3.Position
        ax_feature_bar.OuterPosition = sp3.OuterPosition

        ax_feature_scatter.Position      = sp4.Position
        ax_feature_scatter.OuterPosition = sp4.OuterPosition
        
        clf(fh_output)
        
        copyobj([fh_logrank.Children],fh_output)
        copyobj([fh_tsne.Children],fh_output)
        copyobj([fh_feature_bar.Children],fh_output)
        copyobj([fh_feature_scatter.Children],fh_output)
     

        close([fh_feature_bar,fh_feature_scatter,fh_tsne,fh_logrank]);
        
        
end
end
