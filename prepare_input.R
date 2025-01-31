library(dplyr)
library(tidyverse)
library(scales)

columns <- list('Diversity'=c('nclusters_final_downsampled',
                              'nclusters_eff_downsampled',
                              'mean_cluster_size_downsampled',
                              'std_cluster_size_downsampled',
                              'cluster_gini_downsampled',
                              'cluster_entropy_downsampled',
                              'cluster_inverse_simpson_downsampled',
                              'top_clone_occupancy_downsampled',
                              'big_clones_occupancy_downsampled'),
                'Isotypes'=c('alpha_ratio_clusters',
                             'frac_clusters_SA1','frac_clusters_SA2','frac_clusters_SE',
                             'frac_clusters_SG1','frac_clusters_SG2','frac_clusters_SG3',
                             'frac_clusters_SG4','frac_clusters_SM',
                             'frac_clusters_Sa','frac_clusters_Sg3',
                             'frac_clusters_Sg2b','frac_clusters_Sg2c',
                             'frac_clusters_Sm','frac_clusters_Sg1'),
                'Structural'=c('insert_frequency',
                               'frac_breaks_inversions','frac_breaks_duplications',
                               'median_duplication_size','median_inversion_size'),
                'Context'=c('homology_fw','homology_rv',
                            'donor_score_S', 'donor_score_W', 
                            'donor_score_WGCW', 'donor_score_GAGCT',
                            'donor_score_GGGST',
                            'donor_complexity',
                            'acceptor_score_S', 'acceptor_score_W', 
                            'acceptor_score_WGCW', 'acceptor_score_GAGCT',
                            'acceptor_score_GGGST',
                            'acceptor_complexity',
                            'n_homology_switch', 
                            'n_untemplated_switch',
                            'frac_blunt'),
                'Breaks'=c('mean_length',
                           'std_length',
                           'mean_GC',
                           'std_GC',
                           'mean_length_SA',
                           'mean_length_SG',
                           'std_length_SA',
                           'std_length_SG',
                           'frac_breaks_single',
                           'frac_breaks_sequential',
                           'frac_breaks_intra',
                           'mean_intra_break_size',
                           'std_intra_break_size',
                           'spread_SM','spread_SA','spread_SG',
                           'spread_SA1','spread_SA2',
                           'spread_SG1','spread_SG2',
                           'spread_SG2B','spread_SG2C',
                           'spread_SG3','spread_SG4',
                           'frac_breaks_SM_SM','frac_breaks_SM_SG',
                           'frac_breaks_SG_SG','frac_breaks_SM_SA',
                           'frac_breaks_SG_SA','frac_breaks_SA_SA'),
                'Context_isotype'=c('homology_fw_SM_SA','homology_rv_SM_SA',
                                    'homology_fw_SM_SG','homology_rv_SM_SG',
                                    'donor_score_GAGCT_SM_SG','donor_score_GGGST_SM_SA',
                                    'donor_score_GGGST_SM_SG','donor_score_S_SM_SA',
                                    'donor_score_S_SM_SG','donor_score_W_SM_SA',
                                    'donor_score_W_SM_SG','donor_score_WGCW_SM_SA',
                                    'donor_score_WGCW_SM_SG','donor_score_GAGCT_SM_SA',
                                    'donor_complexity_SM_SG','donor_complexity_SM_SA',
                                    'acceptor_score_GAGCT_SM_SA','acceptor_score_GAGCT_SM_SG',
                                    'acceptor_score_GGGST_SM_SA','acceptor_score_GGGST_SM_SG',
                                    'acceptor_score_S_SM_SA','acceptor_score_WGCW_SM_SG',  
                                    'acceptor_score_S_SM_SG','acceptor_score_W_SM_SA',
                                    'acceptor_score_W_SM_SG','acceptor_score_WGCW_SM_SA',
                                    'acceptor_complexity_SM_SA',
                                    'acceptor_complexity_SM_SG',
                                    'n_homology_switch_SA', 
                                    'n_homology_switch_SG', 
                                    'n_untemplated_switch_SA',
                                    'n_untemplated_switch_SG',
                                    'frac_blunt_SA',
                                    'frac_blunt_SG'),
                'Diversity_raw'=c('nclusters_final',
                                  'nclusters_eff',
                                  'mean_cluster_size',
                                  'std_cluster_size',
                                  'cluster_gini',
                                  'cluster_entropy',
                                  'cluster_inverse_simpson',
                                  'top_clone_occupancy',
                                  'big_clones_occupancy'),
                'Variants'=c('frac_variants_germline',
                             'frac_variants_transitions',
                             'num_variants')) %>%
  unlist() %>%
  data.frame(col=.) %>%
  tibble::rownames_to_column('category') %>%
  dplyr::mutate(category=factor(gsub('[0-9]*$','',category), 
                                levels=c('Isotypes',
                                         'Diversity',
                                         'Diversity_raw',
                                         'Structural',
                                         'Context',
                                         'Context_isotype',
                                         'Breaks',
                                         'Variants')))

# synthetic data
synthetic.data <- read.csv('../benchmarks/20240313_cluster_benchmark_dense/output/summary/20240313_cluster_benchmark_dense_stats.csv',header=1,row.names=1) %>%
  tibble::rownames_to_column('sample') %>%
  dplyr::filter(nreads_final >= 500) %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(nclones=as.numeric(gsub('mix_n([0-9]*)_s([0-9]*)_([0-9])*k','\\1',sample)),
                rep=gsub('mix_n([0-9]*)_s([0-9]*)_([0-9])*k','\\2',sample),
                nreads=gsub('mix_n([0-9]*)_s([0-9]*)_([0-9]*k)','\\3',sample)) %>%
  dplyr::select(tidyselect::any_of(c('sample','nclones','rep','nreads',columns$col))) 
write.csv(synthetic.data, 'data/synthetic_data.csv')

# human and mouse results

human_results <- read.csv('../human_samples/20240307_hg38_results/output/summary/20240307_human_samples_stats.csv',
                          header=1,row.names=1)
human_results_by_reads <- read.csv('../human_samples/20240923_reads_averaging/output/summary/20240923_human_samples_stats.csv',
                                   header=1,row.names=1)

mouse_results <- read.csv('../mouse_samples/20240314_results/output/summary/20240314_mouse_samples_stats.csv',
                          header=1,row.names=1) 
mouse_results_no_Sg <- read.csv('../mouse_samples/20240905_no_SG/output/summary/20240905_mouse_samples_stats.csv',
                                header=1,row.names=1)

# cell switch
cellswitch_sample_sheet <- read.csv('../sodar/2020_CellSwitch/a_2020_CellSwitch.txt',
                            sep='\t',header=1) %>%
  dplyr::left_join(read.csv('../sodar/2020_CellSwitch/s_2020_CellSwitch.txt',sep='\t',header=1),
                   by='Sample.Name') %>%
  dplyr::mutate(Donor=as.character(Source.Name),
                sample=Sample.Name,
                ncells=Parameter.Value.Cell.Number.,
                Batch=as.character(Parameter.Value.Batch.),
                QC=as.character(Parameter.Value.QC.)) %>%
  dplyr::select(Donor,sample,ncells,Batch,QC) %>%
  dplyr::filter(Donor %in% c('21084','21085','21086'))

cellswitch.data <- human_results %>%
  tibble::rownames_to_column('sample') %>%
  dplyr::inner_join(cellswitch_sample_sheet, by='sample') %>%
  dplyr::filter(nreads_final >= 500, QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(Donor=plyr::revalue(Donor, c('21084'='A',
                                             '21085'='B',
                                             '21086'='C'))) %>%
  dplyr::select(Donor,ncells,Batch,sample,nclusters_final,alpha_ratio_clusters)
write.csv(cellswitch.data, 'data/cellswitch_data.csv')

# HD Cohort 1
HD_sample_sheet <- rbind(read.csv('../sodar/2023_OC/a_2023_OC.txt',
                                  sep='\t',header=1) %>%
                           dplyr::left_join(read.csv('../sodar/2023_OC/s_2023_OC.txt',sep='\t',header=1),
                                            by='Sample.Name') %>%
                           dplyr::mutate(Donor=as.character(Source.Name),
                                         Sample=Extract.Name.1,
                                         Batch=as.character(Parameter.Value.Batch.),
                                         QC=Parameter.Value.QC.,
                                         Age=Characteristics.Age.,
                                         Sex=Characteristics.Sex.,
                                         Cohort='C1',
                                         Sequencing=Parameter.Value.Sequencing.method.) %>%
                           dplyr::filter(grepl('HD',Donor)) %>%
                           dplyr::select(Donor,Sample,Batch,Sequencing,QC,Age,Sex,Cohort),
                         read.csv('../sodar/2023_HD/a_2023_HD.txt',
                                  sep='\t',header=1) %>%
                           dplyr::left_join(read.csv('../sodar/2023_HD/s_2023_HD.txt',sep='\t',header=1),
                                            by='Sample.Name') %>%
                           dplyr::mutate(Donor=Source.Name,
                                         Batch=as.character(Parameter.Value.Batch.),
                                         QC=Parameter.Value.QC.,
                                         Sample=Extract.Name.1,
                                         Age=Characteristics.Age.,
                                         Sex=Characteristics.Sex.,
                                         Cohort='C1',
                                         Sequencing=Parameter.Value.Sequencing.method.) %>%
                           dplyr::filter(grepl('HD',Donor)) %>%
                           dplyr::select(Donor,Batch,Sample,Sequencing,QC,Age,Sex,Cohort),
                         read.csv('../sodar/2023_MS/a_2023_MS.txt',
                                  sep='\t',header=1) %>%
                           dplyr::left_join(read.csv('../sodar/2023_MS/s_2023_MS.txt',sep='\t',header=1),
                                            by='Sample.Name') %>%
                           dplyr::mutate(Donor=Source.Name,
                                         Batch=as.character(Parameter.Value.Batch.),
                                         QC=Parameter.Value.QC.,
                                         Age=Characteristics.Age.,
                                         Sex=Characteristics.Sex.,
                                         Sample=Extract.Name.1,
                                         Cohort='C1',
                                         Sequencing=Parameter.Value.Sequencing.method.) %>%
                           dplyr::filter(grepl('HD',Donor)) %>%
                           dplyr::select(Donor,Batch,Sample,Sequencing,QC,Age,Sex,Cohort),
                           read.csv('../sodar/2023_FB/a_2023_FB.txt',
                                  sep='\t',header=1) %>%
                           dplyr::left_join(read.csv('../sodar/2023_FB/s_2023_FB.txt',sep='\t',header=1),
                                            by='Sample.Name') %>%
                           dplyr::mutate(Donor=Source.Name,
                                         Batch=as.character(Parameter.Value.Batch.),
                                         QC=Parameter.Value.QC.,
                                         Age=Parameter.Value.Age.,
                                         Sex=Characteristics.Sex.,
                                         Diagnosis=Characteristics.Diagnosis.,
                                         Sample=Extract.Name.1,
                                         Cohort='C2',
                                         Sequencing=Parameter.Value.Sequencing.method.) %>%
                           dplyr::filter(Diagnosis=='control',grepl('2023',Sample)) %>%
                           dplyr::select(Donor,Batch,Sample,Sequencing,QC,Age,Sex,Cohort))

HD.donors <- unique(HD_sample_sheet$Donor[HD_sample_sheet$Cohort=='C1'])
new.donors <- paste0('donor',seq(1:length(HD.donors)))

HD.samples <- unique(HD_sample_sheet$Sample[HD_sample_sheet$Cohort=='C1'])
new.samples <- paste0('sample',seq(1:length(HD.samples)))
feature_map <- setNames(seq(1:length(columns$col)),columns$col)

HD.data <- human_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::inner_join(HD_sample_sheet, by='Sample') %>%
  dplyr::filter(nreads_final >= 500, QC=='PASS') %>%
  dplyr::select(tidyselect::any_of(c('Donor','Batch','Sample','Sequencing','Age','Sex','Cohort',columns$col))) %>%
  dplyr::mutate(Age=ifelse(Cohort=='C1',ifelse(Age < 50,'<50','>=50'),Age)) %>%
  gather(feature,value,-c(Donor,Batch,Sample,Sequencing,Age,Sex,Cohort)) %>%
  dplyr::group_by(feature,Cohort) %>%
  dplyr::mutate(Donor=ifelse(Cohort=='C1',
                             paste0(plyr::revalue(Donor, setNames(sample(new.donors,n_distinct(Donor)),unique(Donor))),'_',feature_map[feature]),
                             Donor),
                Sample=ifelse(Cohort=='C1',
                              paste0(plyr::revalue(Sample, setNames(sample(new.samples,n_distinct(Sample)),unique(Sample))),'_',feature_map[feature]),
                              Sample)) %>%
  dplyr::ungroup()
write.csv(HD.data, 'data/HD_data.csv')

HD.data.by_reads <- human_results_by_reads %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::inner_join(HD_sample_sheet, by='Sample') %>%
  dplyr::filter(nreads_final >= 500, QC=='PASS') %>%
  dplyr::select(tidyselect::any_of(c('Donor','Batch','Sample','Sequencing','Age','Sex','Cohort',columns$col))) %>%
  dplyr::mutate(Age=ifelse(Cohort=='C1',ifelse(Age < 50,'<50','>=50'),Age)) %>%
  gather(feature,value,-c(Donor,Batch,Sample,Sequencing,Age,Sex,Cohort)) %>%
  dplyr::group_by(feature,Cohort) %>%
  dplyr::mutate(Donor=ifelse(Cohort=='C1',
                             paste0(plyr::revalue(Donor, setNames(sample(new.donors,n_distinct(Donor)),unique(Donor))),'_',feature_map[feature]),
                             Donor),
                Sample=ifelse(Cohort=='C1',
                              paste0(plyr::revalue(Sample, setNames(sample(new.samples,n_distinct(Sample)),unique(Sample))),'_',feature_map[feature]),
                              Sample))
write.csv(HD.data.by_reads, 'data/HD_data_by_reads.csv')

CH12_sample_sheet <- read.csv('../sodar/2022_CH12/a_2022_CH12.txt',
                              sep='\t',header=1) %>%
  dplyr::left_join(read.csv('../sodar/2022_CH12/s_2022_CH12.txt',sep='\t',header=1),
                   by='Sample.Name') %>%
  dplyr::mutate(Genotype=Characteristics.Genotype.,
                Batch=as.character(Parameter.Value.Batch.),
                QC=Parameter.Value.QC.,
                Sample=Extract.Name.1) %>%
  dplyr::select(Genotype,Batch,Sample,QC)

CH12.data <- mouse_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(CH12_sample_sheet, by='Sample') %>%
  dplyr::filter(nreads_final >= 500,!grepl('PP|exon|noAct|20220411|20220810|20220618',Sample),QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Genotype','Batch','Sample',columns$col))) 
write.csv(CH12.data, 'data/CH12_data.csv')

CH12.data.noSg <- mouse_results_no_Sg %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(CH12_sample_sheet, by='Sample') %>%
  dplyr::filter(nreads_final >= 500,!grepl('PP|exon|noAct|20220411|20220810|20220618',Sample),QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Genotype','Batch','Sample',columns$col))) 
write.csv(CH12.data.noSg, 'data/CH12_data_noSg.csv')

mouse_sample_sheet <- read.csv('../sodar/2019_tests/a_2019_test.txt',
                               sep='\t',header=1) %>%
  dplyr::left_join(read.csv('../sodar/2019_tests/s_2019_test.txt',sep='\t',header=1),
                   by='Sample.Name') %>%
  dplyr::mutate(material=gsub(" KO",'',Parameter.Value.Template.),
                Batch=as.character(Parameter.Value.Batch.),
                QC=Parameter.Value.QC.,
                Sample=Extract.Name.1,
                Genotype='WT') %>%
  dplyr::filter(Batch=='20240817') %>%
  dplyr::select(material,Batch,Genotype,Sample,QC)
  
                            
mouse.data <- mouse_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(mouse_sample_sheet, by='Sample') %>%
  dplyr::filter(nreads_final >= 500) %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('material','Batch','Sample',columns$col)))
write.csv(mouse.data, 'data/mouse_data.csv')

mouse.data.noSg <- mouse_results_no_Sg %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(mouse_sample_sheet, by='Sample') %>%
  dplyr::filter(nreads_final >= 500) %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('material','Batch','Sample',columns$col))) 
write.csv(mouse.data.noSg, 'data/mouse_data_noSg.csv')

QP_sample_sheet <- read.csv('../sodar/2019_QianPan/a_2019_QianPan.txt',
                            sep='\t',header=1) %>%
  dplyr::left_join(read.csv('../sodar/2019_QianPan/s_2019_QianPan.txt',sep='\t',header=1),
                   by='Sample.Name') %>%
  dplyr::mutate(Donor=Source.Name,
                Diagnosis=Characteristics.Genotype.,
                Sex=Characteristics.Sex.,
                Age=Characteristics.Age.,
                Batch=as.character(Parameter.Value.Batch.),
                QC=Parameter.Value.QC.,
                Sample=Extract.Name.1) %>%
  dplyr::select(Donor,Diagnosis,Sex,Age,Batch,QC,Sample)

QP.data <- human_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(QP_sample_sheet, by='Sample') %>%
  dplyr::filter(nreads_final > 500, QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Donor','Batch','Sample','Diagnosis','Age','Sex',columns$col)))
write.csv(QP.data, 'data/QP_data.csv')

FB_sample_sheet <- read.csv('../sodar/2023_FB/a_2023_FB.txt',
                            sep='\t',header=1) %>%
  dplyr::left_join(read.csv('../sodar/2023_FB/s_2023_FB.txt',sep='\t',header=1),
                   by='Sample.Name') %>%
  dplyr::mutate(Donor=Source.Name,
                Sample_donor=Sample.Name,
                Diagnosis=plyr::revalue(ifelse(Donor %in% c('FB_13','FB_21'),'Atypical',Characteristics.Diagnosis.),
                                        c('control'='Control',
                                          'subclass'='Subclass')),
                complication=Characteristics.Complication.,
                Sex=Characteristics.Sex.,
                splenomegaly=Characteristics.splenomegaly.,
                lymphoproliferation=Characteristics.lymphoproliferation.,
                lung=Characteristics.lung.,
                PnPS=Characteristics.PnPS.,
                Tetanus=Characteristics.Tetanus.,
                Diphteria=Characteristics.Diphteria.,
                Timepoint=Parameter.Value.timepoint.,
                Age=Parameter.Value.Age.,
                EUROclass=Parameter.Value.EUROclass.,
                pct_B=Parameter.Value.pct_B.,
                pct_IgA=Parameter.Value.pct_IgA.,
                pct_IgG=Parameter.Value.pct_IgG.,
                Batch=as.character(Parameter.Value.Batch.),
                QC=Parameter.Value.QC.,
                Sample=Extract.Name.1) %>%
  dplyr::distinct(Donor,Sample_donor,Diagnosis,QC,Sample,.keep_all = TRUE)

FB.data <- human_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(FB_sample_sheet %>% dplyr::distinct(Donor, Diagnosis, Timepoint, Sample_donor, Sex, Age, Sample, QC, Batch), by='Sample') %>%
  dplyr::filter(nreads_final >= 500, QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Donor','Diagnosis','Timepoint','Sample_donor','Sex','Age','Sample','Batch',columns$col))) 
write.csv(FB.data, 'data/FB_data.csv')

FB.data.pooled <- human_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::filter(grepl('pooled_',Sample)) %>%
  dplyr::mutate(Sample_donor=gsub('pooled_','',Sample)) %>%
  dplyr::inner_join(FB_sample_sheet %>% dplyr::distinct(Donor, Diagnosis, Timepoint, Sample_donor, Sex, Age, 
                                                        splenomegaly, lymphoproliferation, Tetanus, Diphteria, PnPS, lung,
                                                        complication, EUROclass, pct_B, pct_IgA, pct_IgG), by='Sample_donor') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(num_variants_somatic=num_variants * (1 - frac_variants_germline)) %>%
  dplyr::select(tidyselect::any_of(c('Donor','Diagnosis','Timepoint','Sample_donor','Sex','Age','Sample',
                                     'Tetanus','Diphteria','PnPS','lung',
                                     'pct_B','pct_IgA','pct_IgG',columns$col)))
write.csv(FB.data.pooled, 'data/FB_data_pooled.csv')

vincendeau.data <- read.csv('../mouse_samples/20240425_vincendeau/output/summary/20240425_vincendeau_samples_stats.csv',
                            header=1,row.names=1) %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(Genotype=gsub('PRJNA[0-9]*_(.*)_[0-9]$','\\1',Sample)) %>%
  dplyr::select(tidyselect::any_of(c('Genotype','Sample',columns$col))) 
write.csv(vincendeau.data, 'data/vincendeau_data.csv')

HTGTS.samples <- c("SRR2104731","SRR2104732","SRR2104733","SRR2104734","SRR2104735","SRR2104736","SRR2104737","SRR2104738",
                   "SRR2104739","SRR2104740","SRR2104741","SRR2104742","SRR2104743","SRR2104744",
                   "SRR6293456","SRR6293457","SRR6293458","SRR6293459","SRR6293460","SRR6293461","SRR6293462","SRR6293463",
                   "SRR6293464","SRR6293465","SRR6293466","SRR6293467","SRR6293468","SRR6293469","SRR6293470","SRR6293471",
                   "SRR6293472","SRR6293473","SRR6293474","SRR6293475","SRR6293476","SRR6293477","SRR6293478","SRR6293479")

HTGTS.meta <- data.frame(celltype=c(rep('spleen',14),rep('CH12',24)),
                         Genotype=c(rep('WT',5),rep('Atm',3),rep('Trp53bp1',6),
                                    rep('WT',6),rep('Atm',6),rep('Trp53bp1',6),rep('Lig4',6)),
                         Sample=HTGTS.samples) %>%
  dplyr::mutate(Genotype=factor(Genotype, levels=c('WT','Atm','Trp53bp1','Lig4')))

HTGTS.data <- read.csv('../mouse_samples/20241204_HTGTS/output/summary/HTGTS_stats.csv',
                       header=1,row.names=1) %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(HTGTS.meta, by='Sample') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Genotype','celltype','Sample',columns$col))) 
write.csv(HTGTS.data, 'data/HTGTS_data.csv')

HTGTS.homology <- lapply(HTGTS.samples, function(smp) read.csv(paste0('../mouse_samples/20241204_HTGTS/pipeline/',smp,'/',smp,'_cluster_analysis.csv'), header=1, row.names=1) %>%
                           tibble::rownames_to_column('cluster') %>%
                           dplyr::filter(cluster %in% (read.csv(paste0('../mouse_samples/20241204_HTGTS/pipeline/',smp,'/',smp,'_clustering.csv'),header=1,row.names=1) %>%
                                                         dplyr::filter(filtered_cluster >= 0)  %>%
                                                         dplyr::pull(unique(filtered_cluster)))) %>%
                           dplyr::select(isotype,n_homology_switch) %>%
                           dplyr::mutate(n_homology=cut(n_homology_switch, breaks=seq(-1,10),
                                                        labels=seq(0,10))) %>%
                           dplyr::group_by(isotype, n_homology) %>%
                           dplyr::summarise(n=n()) %>%
                           dplyr::group_by(isotype) %>%
                           dplyr::mutate(frac=n/sum(n),
                                         Sample=smp)) %>%
  do.call(rbind,.) %>%
  dplyr::inner_join(HTGTS.meta, by='Sample') %>%
  dplyr::filter(!is.na(n_homology))
write.csv(HTGTS.homology, 'data/HTGTS_homology.csv')

vdj_samples <- read.csv('../human_samples/20230105_VDJ/20230113_mixcr/samples.txt',sep='\t',header=FALSE)[,2]
VDJ_results <- lapply(vdj_samples, function(sample)
  read.csv(file.path('..','human_samples','20230105_VDJ','20230113_mixcr',
                     sample,paste0(sample,'.clones_IGH.tsv')),
           sep='\t',header=1) %>%
    dplyr::mutate(Sample=sample)) %>%
  do.call(rbind,.) %>%
  dplyr::mutate(Donor=plyr::revalue(as.character(gsub('_.*','',Sample)),
                                    c('21084'='A',
                                      '21085'='B',
                                      '21086'='C')),
                ncells=as.numeric(gsub('^[0-9]*_([0-9]*)_([0-9])','\\1',Sample))) %>%
  dplyr::select(cloneId,readCount,uniqueMoleculeCount,allVHitsWithScore,allCHitsWithScore,aaSeqCDR3,Sample,Donor,ncells)
write.csv(VDJ_results, 'data/VDJ_results.csv')

ont <- read.csv("../mouse_samples/demultiplexing/barcodes_primers.fa",header=FALSE, sep='\t')
barcodes <- setNames(ont[!grepl('^>',ont$V1),'V1'],gsub('^>','',ont[grepl("^>",ont$V1),'V1']))

GEO_sample_sheet <- rbind(read.csv('../sodar/2022_CH12/a_2022_CH12.txt',
                                   sep='\t',header=1) %>%
                            dplyr::left_join(read.csv('../sodar/2022_CH12/s_2022_CH12.txt',sep='\t',header=1),
                                             by='Sample.Name') %>%
                            dplyr::mutate(Genotype=Characteristics.Genotype.,
                                          Batch=as.character(Parameter.Value.Batch.),
                                          material='CH12',
                                          barcode=Parameter.Value.Barcode.,
                                          PCR_cycles=Parameter.Value.PCR.cycles.,
                                          PCR_Volume=Parameter.Value.PCR.volume.,
                                          primers=Parameter.Value.Primers.,
                                          barcode=Parameter.Value.Barcode.,
                                          Polymerase=Parameter.Value.DNA.polymerase.,
                                          Activation=Parameter.Value.activation.,
                                          Sample=Extract.Name.1) %>%
                            dplyr::filter(!grepl("20220411|20220810|20220618",Batch)) %>%
                            dplyr::select(Sample,material,Genotype,barcode,PCR_cycles,PCR_Volume,primers,barcode,Polymerase,Activation,Batch),
                          read.csv('../sodar/2019_tests/a_2019_test.txt',
                                   sep='\t',header=1) %>%
                            dplyr::left_join(read.csv('../sodar/2019_tests/s_2019_test.txt',sep='\t',header=1),
                                             by='Sample.Name') %>%
                            dplyr::mutate(material=gsub(" KO",'',Parameter.Value.Template.),
                                          Batch=as.character(Parameter.Value.Batch.),
                                          PCR_cycles=Parameter.Value.PCR.cycles.,
                                          PCR_Volume=Parameter.Value.PCR.volume.,
                                          primers="M-G,M-A",
                                          barcode=Parameter.Value.Barcode.,
                                          Polymerase=Parameter.Value.DNA.polymerase.,
                                          Activation="",
                                          Sample=Extract.Name.1,
                                          Genotype='WT') %>%
                            dplyr::filter(Batch=='20240817') %>%
                            dplyr::select(Sample,material,Genotype,barcode,PCR_cycles,PCR_Volume,primers,barcode,Polymerase,Activation,Batch)) %>%
  dplyr::group_by(material, Genotype, primers) %>% 
  dplyr::mutate(replicate=seq(1:n()),
                sequence=barcodes[barcode]) %>%
  write.csv('GEO_sample_sheet.csv')
                          
                          