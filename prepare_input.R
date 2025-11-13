library(dplyr)
library(tidyverse)
library(scales)
library(openxlsx)

set.seed(1)

columns <- list('Diversity'=c('clusters',
                              'clusters_eff',
                              'cluster_size_mean',
                              'cluster_size_std',
                              'cluster_gini',
                              'cluster_entropy',
                              'cluster_inverse_simpson',
                              'occupancy_top_cluster',
                              'occupancy_big_clusters'),
                'Isotypes'=c('alpha_ratio_clusters',
                             'pct_clusters_SA1','pct_clusters_SA2','pct_clusters_SE',
                             'pct_clusters_SG1','pct_clusters_SG2','pct_clusters_SG3',
                             'pct_clusters_SG4','pct_clusters_SM',
                             'pct_clusters_Sa','pct_clusters_Sg3',
                             'pct_clusters_Sg2b','pct_clusters_Sg2c',
                             'pct_clusters_Sm','pct_clusters_Sg1'),
                'Structural'=c('pct_templated_inserts',
                               'pct_inversions','pct_duplications',
                               'duplication_size','inversion_size'),
                'Context'=c('homology_score_fw','homology_score_rv',
                            'donor_score_S', 'donor_score_W', 
                            'donor_score_WGCW', 'donor_score_GAGCT',
                            'donor_score_GGGST',
                            'donor_complexity',
                            'acceptor_score_S', 'acceptor_score_W', 
                            'acceptor_score_WGCW', 'acceptor_score_GAGCT',
                            'acceptor_score_GGGST',
                            'acceptor_complexity',
                            'homology', 
                            'untemplated_inserts',
                            'pct_blunt'),
                'Breaks'=c('read_length_mean',
                           'read_length_std',
                           'GC_content_mean',
                           'GC_content_std',
                           'read_length_SA_mean',
                           'read_length_SG_mean',
                           'read_length_SA_std',
                           'read_length_SG_std',
                           'pct_direct_switch',
                           'pct_sequential_switch',
                           'pct_intraswitch_deletion',
                           'intraswitch_size_mean',
                           'intraswitch_size_std',
                           'break_dispersion_SM','break_dispersion_SA','break_dispersion_SG',
                           'break_dispersion_SA1','break_dispersion_SA2',
                           'break_dispersion_SG1','break_dispersion_SG2',
                           'break_dispersion_SG2B','break_dispersion_SG2C',
                           'break_dispersion_SG3','break_dispersion_SG4',
                           'pct_breaks_SM_SM','pct_breaks_SM_SG',
                           'pct_breaks_SG_SG','pct_breaks_SM_SA',
                           'pct_breaks_SG_SA','pct_breaks_SA_SA'),
                'Context_isotype'=c('homology_score_fw_SM_SA','homology_score_rv_SM_SA',
                                    'homology_score_fw_SM_SG','homology_score_rv_SM_SG',
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
                                    'homology_SA', 
                                    'homology_SG', 
                                    'untemplated_inserts_SA',
                                    'untemplated_inserts_SG',
                                    'pct_blunt_SA',
                                    'pct_blunt_SG'),
                'Diversity_raw'=c('clusters_raw',
                                  'clusters_eff_raw',
                                  'cluster_size_mean_raw',
                                  'cluster_size_std_raw',
                                  'cluster_gini_raw',
                                  'cluster_entropy_raw',
                                  'cluster_inverse_simpson_raw',
                                  'occupancy_top_cluster_raw',
                                  'occupancy_big_clusters_raw'),
                'Variants'=c('somatic_variants',
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

wb <- createWorkbook()
addWorksheet(wb, sheetName='description')
addWorksheet(wb, sheetName='synthetic data')
addWorksheet(wb, sheetName='cell numbers')
addWorksheet(wb, sheetName='meta-clustering')
addWorksheet(wb, sheetName='VDJ data')
addWorksheet(wb, sheetName='HD data')
addWorksheet(wb, sheetName='HD data by reads')
addWorksheet(wb, sheetName='Vincendeau data')
addWorksheet(wb, sheetName='HTGTS data')
addWorksheet(wb, sheetName='HTGTS homology')
addWorksheet(wb, sheetName='mouse tissues')
addWorksheet(wb, sheetName='mouse tissues (no Sg)')
addWorksheet(wb, sheetName='CH12 data')
addWorksheet(wb, sheetName='CH12 data (no Sg)')
addWorksheet(wb, sheetName='CVID data')
addWorksheet(wb, sheetName='CVID data pooled')
addWorksheet(wb, sheetName='DNA repair data')

description <- data.frame(
  rbind(c('synthetic data','data from simulation experiments','1D'),
        c('cell numbers','data from runs with different cell numbers','1E, 1F, 1G'),
        c('meta-clustering','meta-clustering of cell number data','1F'),
        c('VDJ data','data from bulk VDJ sequencing','1G'),
        c('HD data','data for cohort C1+C2','2B, 2C, 2D, S2F, S2G, S2I'),
        c('HD data by reads','data for cohort C1+C2 (no clustering)','S2I'),
        c('CVID data','data for CVID cohort + controls','2D, 4B, S4A'),
        c('Vincendeau data','data from Vincendeau et al.','2F, 2G, 2H, S2J'),
        c('HTGTS data','re-processed HTGTS data','S3A-C'),
        c('HTGTS homology','homology results for HTGTS data','S3D'),
        c('mouse tissues','data from mouse lymph nodes + spleen','2F, 3B, S2J'),
        c('mouse tissues (no Sg)','data from mouse lymph nodes + spleen (without Sg reads)','3B'),
        c('CH12 data','results for CH12 WT and KO cells','3B-G'),
        c('CH12 data (no Sg)','CH12 results without Sg reads','3B'),
        c('CVID (pooled)','CVID cohort (pooled replicates)','4C-J, S4B-H'),
        c('DNA repair','DNA repair cohort', '4F, 4G, 4J')))
colnames(description) <- c('sheet','description','used in figure')

writeData(wb, 'description', description)

# synthetic data
synthetic.data <- read.csv('swibrid_runs/benchmark/dense/output/summary/20240313_cluster_benchmark_dense_stats.csv',header=1,row.names=1) %>%
  tibble::rownames_to_column('sample') %>%
  dplyr::filter(reads >= 500) %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(nclones=as.numeric(gsub('mix_n([0-9]*)_s([0-9]*)_([0-9])*k','\\1',sample)),
                rep=gsub('mix_n([0-9]*)_s([0-9]*)_([0-9])*k','\\2',sample),
                nreads=gsub('mix_n([0-9]*)_s([0-9]*)_([0-9]*k)','\\3',sample)) %>%
  dplyr::select(tidyselect::any_of(c('sample','nclones','rep','nreads',columns$col))) 
writeData(wb, 'synthetic data', synthetic.data)

# human and mouse results

human_results <- read.csv('swibrid_runs/human/output/summary/20240307_human_samples_stats.csv',
                          header=1,row.names=1)
human_results_by_reads <- read.csv('swibrid_runs/human/output_read_averaging/summary/20240923_human_samples_stats.csv',
                                   header=1,row.names=1)

mouse_results <- read.csv('swibrid_runs/mouse/output/summary/20240314_mouse_samples_stats.csv',
                          header=1,row.names=1) 
mouse_results_no_Sg <- read.csv('swibrid_runs/mouse/output_no_Sg/summary/20240905_mouse_samples_stats.csv',
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
  dplyr::filter(reads >= 500, QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(Donor=plyr::revalue(Donor, c('21084'='A',
                                             '21085'='B',
                                             '21086'='C'))) %>%
  dplyr::select(Donor,ncells,Batch,sample,clusters_raw,alpha_ratio_clusters,cluster_gini_raw,cluster_inverse_simpson_raw)
writeData(wb, 'cell numbers', cellswitch.data)

meta_clustering <- lapply(c('21084','21085','21086','mixed'), 
                          function(donor) 
                            lapply(if (donor=='mixed') seq(1,5) else seq(1,3),
                                   function(rep)
                                     lapply(c(50000,100000), function(nc) 
                                     {
                                       fname <- paste0('swibrid_runs/human/meta_clustering/',
                                                       donor,'_',rep,'_',format(nc,scientific=FALSE),'_meta_clustering.csv')
                                       if (file.exists(fname)) 
                                         read.csv(fname, header=1,row.names=1) %>%
                                         dplyr::group_by(sample, meta_cluster) %>%
                                         dplyr::summarise(nreads=n(),
                                                          n_clusters=n_distinct(cluster)) %>%
                                         dplyr::ungroup() %>%
                                         dplyr::group_by(meta_cluster) %>%
                                         dplyr::mutate(n_samples=n_distinct(sample[nreads > 0 * sum(nreads)])) %>%
                                         dplyr::ungroup() %>%
                                         dplyr::group_by(sample) %>%
                                         dplyr::mutate(cluster_size=cut(nreads/sum(nreads), c(0,.0001,.001,.01,.1,1), 
                                                                        labels=c('<.01%','<.1%','<1%','<10%','>10%'))) %>%
                                         dplyr::mutate(Donor=donor,
                                                       ncells=nc,
                                                       replicate=rep)
                                     }) %>%
                                     do.call(rbind, .)) %>%
                            do.call(rbind, .)) %>%
  do.call(rbind, .) %>%
  dplyr::group_by(Donor,replicate,ncells) %>%
  dplyr::summarise(frac_clusters=mean(n_samples > 1),
                   frac_ambiguous=mean(n_clusters > 1)) %>%
  dplyr::ungroup() %>% 
  dplyr::distinct(Donor, ncells, frac_clusters, frac_ambiguous) %>%
  dplyr::mutate(mixing=ifelse(Donor=='mixed','inter-donor','intra-donor'),
                Donor=plyr::revalue(Donor, c('21084'='A',
                                             '21085'='B',
                                             '21086'='C')),
                ncells=factor(ifelse(ncells==50000,'50k','100k'),levels=c('50k','100k')))
writeData(wb, 'meta-clustering', meta_clustering)

# HD Cohort 1
HD_sample_sheet <- rbind(read.csv('../sodar/2023_OC/a_2023_OC.txt',
                                  sep='\t',header=1) %>%
                           dplyr::left_join(read.csv('../sodar/2023_OC/s_2023_OC.txt',sep='\t',header=1),
                                            by='Sample.Name') %>%
                           dplyr::mutate(Donor=as.character(Source.Name),
                                         Sample=Comment.Extract.Name.,
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
                                         Sample=Comment.Extract.Name.,
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
                                         Sample=Comment.Extract.Name.,
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
  dplyr::filter(reads >= 500, QC=='PASS') %>%
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
writeData(wb, 'HD data', HD.data)

HD.data.by_reads <- human_results_by_reads %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::inner_join(HD_sample_sheet, by='Sample') %>%
  dplyr::filter(reads >= 500, QC=='PASS') %>%
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
writeData(wb, 'HD data by reads', HD.data.by_reads)

CH12_sample_sheet <- read.csv('../sodar/2022_CH12/a_2022_CH12.txt',
                              sep='\t',header=1) %>%
  dplyr::left_join(read.csv('../sodar/2022_CH12/s_2022_CH12.txt',sep='\t',header=1),
                   by='Sample.Name') %>%
  dplyr::mutate(Genotype=Characteristics.Genotype.,
                Batch=as.character(Parameter.Value.Batch.),
                QC=Parameter.Value.QC.,
                Sample=Comment.Extract.Name.) %>%
  dplyr::select(Genotype,Batch,Sample,QC)

CH12.data <- mouse_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(CH12_sample_sheet, by='Sample') %>%
  dplyr::filter(reads >= 500,!grepl('PP|exon|noAct|20220411|20220810|20220618',Sample),QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Genotype','Batch','Sample',columns$col))) 
writeData(wb, 'CH12 data', CH12.data)

CH12.data.noSg <- mouse_results_no_Sg %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(CH12_sample_sheet, by='Sample') %>%
  dplyr::filter(reads >= 500,!grepl('PP|exon|noAct|20220411|20220810|20220618',Sample),QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Genotype','Batch','Sample',columns$col))) 
writeData(wb, 'CH12 data (no Sg)', CH12.data.noSg)

mouse_sample_sheet <- read.csv('../sodar/2019_tests/a_2019_test.txt',
                               sep='\t',header=1) %>%
  dplyr::left_join(read.csv('../sodar/2019_tests/s_2019_test.txt',sep='\t',header=1),
                   by='Sample.Name') %>%
  dplyr::mutate(material=gsub(" KO",'',Parameter.Value.Template.),
                Batch=as.character(Parameter.Value.Batch.),
                QC=Parameter.Value.QC.,
                Sample=Comment.Extract.Name.,
                Genotype='WT') %>%
  dplyr::filter(Batch=='20240817') %>%
  dplyr::select(material,Batch,Genotype,Sample,QC)
  
                            
mouse.data <- mouse_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(mouse_sample_sheet, by='Sample') %>%
  dplyr::filter(reads >= 500) %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('material','Batch','Sample',columns$col)))
writeData(wb, 'mouse tissues', mouse.data)

mouse.data.noSg <- mouse_results_no_Sg %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(mouse_sample_sheet, by='Sample') %>%
  dplyr::filter(reads >= 500) %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('material','Batch','Sample',columns$col))) 
writeData(wb, 'mouse tissues (no Sg)', mouse.data.noSg)

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
                Sample=Comment.Extract.Name.) %>%
  dplyr::select(Donor,Diagnosis,Sex,Age,Batch,QC,Sample)

QP.data <- human_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(QP_sample_sheet, by='Sample') %>%
  dplyr::filter(reads > 500, QC=='PASS') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Donor','Batch','Sample','Diagnosis','Age','Sex',columns$col)))
writeData(wb, 'DNA repair data', QP.data)

FB_sample_sheet <- read.csv('../sodar/2023_FB/a_2023_FB.txt',
                            sep='\t',header=1) %>%
  dplyr::left_join(read.csv('../sodar/2023_FB/s_2023_FB.txt',sep='\t',header=1),
                   by='Sample.Name') %>%
  dplyr::mutate(Donor=Source.Name,
                Sample_donor=Sample.Name,
                Diagnosis=plyr::revalue(ifelse(Donor %in% c('FB_13','FB_21'),'Atypical',Characteristics.Diagnosis.),
                                        c('control'='Control',
                                          'subclass'='Subclass')),
                Complication=Characteristics.Complication.,
                Sex=Characteristics.Sex.,
                PnPS=Characteristics.PnPS.,
                Tetanus=Characteristics.Tetanus.,
                Diphteria=Characteristics.Diphteria.,
                Timepoint=Parameter.Value.timepoint.,
                Age=Parameter.Value.Age.,
                pct_B=Parameter.Value.pct_B.,
                pct_IgA=Parameter.Value.pct_IgA.,
                pct_IgG=Parameter.Value.pct_IgG.,
                Treatment=Parameter.Value.Treatment.,
                IV=Parameter.Value.IV.,
                Batch=as.character(Parameter.Value.Batch.),
                QC=Parameter.Value.QC.,
                Sample=Comment.Extract.Name.) %>%
  dplyr::distinct(Donor,Sample_donor,Diagnosis,QC,Sample,Complication,.keep_all = TRUE)

FB.data <- human_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(FB_sample_sheet, by='Sample') %>%
  dplyr::filter(reads >= 500, QC=='PASS', grepl('FB', Donor)) %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Donor','Diagnosis','Timepoint','Sample_donor','Sex','Age','Sample','Batch','Treatment','IV','Complication',columns$col))) %>%
  dplyr::mutate(Treatment=ifelse((Donor=='FB_35') | (Treatment==''),'none',Treatment),
                IV=ifelse(IV=='','none',IV),
                Complication=ifelse(Complication=='', 'na', Complication)) 
writeData(wb, 'CVID data', FB.data)

FB.data.pooled <- human_results %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::filter(grepl('pooled_FB',Sample)) %>%
  dplyr::mutate(Sample_donor=gsub('pooled_','',Sample)) %>%
  dplyr::inner_join(FB_sample_sheet %>% dplyr::distinct(Donor, Diagnosis, Timepoint, Sample_donor, Sex, Age, 
                                                        Tetanus, Diphteria, PnPS, Treatment, IV, Complication,
                                                        pct_B, pct_IgA, pct_IgG), by='Sample_donor') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(num_variants_somatic=num_variants * (1 - frac_variants_germline)) %>%
  dplyr::select(tidyselect::any_of(c('Donor','Diagnosis','Timepoint','Sample_donor','Sex','Age','Sample',
                                     'Tetanus','Diphteria','PnPS','Treatment','IV','Complication',
                                     'pct_B','pct_IgA','pct_IgG',columns$col))) %>%
  dplyr::mutate(Treatment=ifelse((Donor=='FB_35') | (Treatment==''),'none',Treatment),
                IV=ifelse(IV=='','none',IV),
                Complication=ifelse(Complication=='', 'na', Complication))
writeData(wb, 'CVID data pooled', FB.data.pooled)

vincendeau.data <- read.csv('swibrid_runs/external/vincendeau/output/summary/20240425_vincendeau_samples_stats.csv',
                            header=1,row.names=1) %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(Genotype=gsub('PRJNA[0-9]*_(.*)_[0-9]$','\\1',Sample)) %>%
  dplyr::select(tidyselect::any_of(c('Genotype','Sample',columns$col))) 
writeData(wb, 'Vincendeau data', vincendeau.data)

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

HTGTS.data <- read.csv('swibrid_runs/external/HTGTS/output/summary/HTGTS_stats.csv',
                       header=1,row.names=1) %>%
  tibble::rownames_to_column('Sample') %>%
  dplyr::inner_join(HTGTS.meta, by='Sample') %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::select(tidyselect::any_of(c('Genotype','celltype','Sample',columns$col))) 
writeData(wb, 'HTGTS data', HTGTS.data)

HTGTS.homology <- lapply(HTGTS.samples, function(smp) read.csv(paste0('swibrid_runs/external/HTGTS/pipeline/',smp,'/',smp,'_cluster_analysis.csv'), header=1, row.names=1) %>%
                           tibble::rownames_to_column('cluster') %>%
                           dplyr::filter(cluster %in% (read.csv(paste0('swibrid_runs/external/HTGTS/pipeline/',smp,'/',smp,'_clustering.csv'),header=1,row.names=1) %>%
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
writeData(wb, 'HTGTS homology', HTGTS.homology)

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
writeData(wb, 'VDJ data', VDJ_results)

saveWorkbook(wb, file = "data/input_data.xlsx", overwrite = TRUE)
