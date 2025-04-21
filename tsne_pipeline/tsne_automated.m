clear
close all
format long
%% Loading table and getting rid of samples with no transcriptional data [preprocessing]
t =readtable('./20180511_cholbiosyn_no_ggps.csv');

used_vars = t.Properties.VariableNames;
gen_var_start = find( strcmp('samples', used_vars))+1;
gen_var_end   = find( strcmp('x_EVENT', used_vars))-1;
gen_vars = used_vars(gen_var_start:gen_var_end);


t = t(~isnan(t{:,gen_var_start}),:);




%% Dividing leukemia samples into peripheral and bone marrow samples (renaming the leukemia disease variable) [preprocessing]

peripheral =  strcmp(t.sample_type_samples,  'Primary_Blood_Derived_Cancer_-_Peripheral_Blood');
peripheral_string_array = repmat({'peripheral'},size(t,1),1);
t.disease_type(peripheral) = join([t.disease_type(peripheral),peripheral_string_array(peripheral)],'_');




marrow =  strcmp(t.sample_type_samples,  'Primary_Blood_Derived_Cancer_-_Bone_Marrow');
marrow_string_array = repmat({'marrow'},size(t,1),1);
t.disease_type(marrow) = join([t.disease_type(marrow),marrow_string_array(marrow)],'_');




%% Filtering for only samples we want [preprocessing]


t = t(strcmp('Solid_Tissue_Normal',t.sample_type_samples)| ...
    strcmp('Primary_Tumor',t.sample_type_samples)| ...
    strcmp('Primary_Blood_Derived_Cancer_-_Peripheral_Blood',t.sample_type_samples)| ...
    strcmp('Primary_Blood_Derived_Cancer_-_Bone_Marrow',t.sample_type_samples),:);



%% Exponentiating and decrementing by one each sample and normalizing to sum of all genes in sample.  RELEVANT
gen_ra = t{:,gen_var_start:gen_var_end};           %  selecting the variables with genetic data
gen_ra = (2.^gen_ra) -(ones(size(gen_ra)));        %  the data were stored incremented by 1, then log-base-two'd
gen_ra = gen_ra ./ repmat(sum(gen_ra,2),1,size(gen_ra,2)); %normalizing the values in each sample to the sum across the sample

t{:,gen_var_start:gen_var_end} =gen_ra;  %rewriting the normalized values to the table
% %


all_diseases = unique(t.disease_type);  %getting this cell array to possibly loop through


%annotation of what is in the cell so I can test on individual diseases

%     1{'Acute_Myeloid_Leukemia_marrow'                                   }
%     2{'Acute_Myeloid_Leukemia_peripheral'                               }
%     3{'Adrenocortical_Carcinoma'                                        }
%     4{'Bladder_Urothelial_Carcinoma'                                    }
%     5{'Brain_Lower_Grade_Glioma'                                        }
%     6{'Breast_Invasive_Carcinoma'                                       }
%     7{'Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma'}
%     8{'Cholangiocarcinoma'                                              }
%     9{'Colon_Adenocarcinoma'                                            }
%     10{'Esophageal_Carcinoma'                                            }
%     11{'Glioblastoma_Multiforme'                                         }
%     12{'Head_and_Neck_Squamous_Cell_Carcinoma'                           }
%     13{'High-Risk_Wilms_Tumor'                                           }
%     14{'Kidney_Chromophobe'                                              }
%     15{'Kidney_Renal_Clear_Cell_Carcinoma'                               }
%     16{'Kidney_Renal_Papillary_Cell_Carcinoma'                           }
%     17{'Liver_Hepatocellular_Carcinoma'                                  }
%     18{'Lung_Adenocarcinoma'                                             }
%     19{'Lung_Squamous_Cell_Carcinoma'                                    }
%     20{'Lymphoid_Neoplasm_Diffuse_Large_B-cell_Lymphoma'                 }
%     21{'Mesothelioma'                                                    }
%     22{'Neuroblastoma'                                                   }
%     23{'Ovarian_Serous_Cystadenocarcinoma'                               }
%     24{'Pancreatic_Adenocarcinoma'                                       }
%     25{'Pheochromocytoma_and_Paraganglioma'                              }
%     26{'Prostate_Adenocarcinoma'                                         }
%     27{'Rectum_Adenocarcinoma'                                           }
%     28{'Sarcoma'                                                         }
%     29{'Skin_Cutaneous_Melanoma'                                         }
%     30{'Stomach_Adenocarcinoma'                                          }
%     31{'Testicular_Germ_Cell_Tumors'                                     }
%     32{'Thymoma'                                                         }
%     33{'Thyroid_Carcinoma'                                               }
%     34{'Uterine_Carcinosarcoma'                                          }
%     35{'Uterine_Corpus_Endometrial_Carcinoma'                            }
%     36{'Uveal_Melanoma'                                                  }


for d= 17 %17 is Liver Hepatocellular Carcinoma
    
    disease_table = t(strcmp(all_diseases{d},t.disease_type),:);
    disease_table = sortrows(disease_table,{'sample'});
    disease_gen_table = disease_table(:,gen_var_start:gen_var_end);
    
    
    %% Writing input that I used in Tensorboard: metadata in first two columns.  I ran the tsne in tensorboard with 'spherize data' checked, perplexity = 14, learning rate = 10.
    writetable(cell2table([num2cell(double(strcmp(disease_table.sample_type_samples,'Solid_Tissue_Normal'))), ...
        disease_table.sample, ...
        num2cell(disease_gen_table{:,:})]),'./tf_input/mathworks_tsne_test.csv','WriteVariableNames',false)
    
    
    
    
    %% Here is the essential_inputdata
    ra_in = disease_gen_table{:,:};
    
    
    
    %% This section 'spherizes' the data.  Finds centroid, shifts, normalizes so that each point is one unit away from origin in multidimensional space
    % I 'sphereized' the data for the attached tensorflow output as
    % well.   You can comment out line 125 to not spherize
    
    
    centroid_new = mean(ra_in,1);   %finding the centroid
    ra_in = ra_in-repmat(centroid_new,size(ra_in,1),1); %shifting data by centroid
    norm_new = repmat(  sqrt(sum(ra_in.^2,2))  , 1,size(ra_in,2) ); %denominator for unit normalization from origin
    ra_in = ra_in./norm_new;  %unit normalizing
    
    
    
    
    disease_table{:,gen_var_start:gen_var_end} = ra_in; %putting the centered normed array back into the table
    
    
    num_clust = 1:10;  %the number of clusters to consider when deciding how many there are
    
    
    % used for automatically selecting the clusters
    cluster_criterion={
        'CalinskiHarabasz'
        'DaviesBouldin'
        'gap'
        'silhouette'
        };
    
    cluster_criterion = cluster_criterion{4};
    
    
    disp(['Running tSNE for ', all_diseases{d}]);
    
    
    %% Setting up and running tsne
    maxiter =5000;
    opts = statset('OutputFcn',@(optimValues,state) tsne_outfun(optimValues,state,disease_table,num_clust,cluster_criterion),'MaxIter',maxiter,'Display','iter','UseParallel',false);
    
    
    
    y=tsne_mod(ra_in,'Algorithm','exact','Perplexity',12,'LearnRate',10,'Options',opts,'Standardize',false,...
        'Distance','euclidean','Verbose',2,'Options',opts,'NumPrint',100,'NumDimensions',2,'Exaggeration',10);  %only difference betwen version from mathworks is the bug I found and corrected
                                                                                                                %same perplexity as was done in tensorboard
    
    
    
end

