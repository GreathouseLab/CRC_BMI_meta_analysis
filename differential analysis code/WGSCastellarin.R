source('~/Dropbox/Workspace/MayoClinic/Stats.R')
load_package()

load("~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/WGS.wk.RData")

data.obj0 <- Castellarin.data.obj.wk
dist.obj0 <- Castellarin.dist.obj.wk

#BMI <- data.obj0$meta.dat$BMI
#BMI.cat <- factor((BMI > 25) + (BMI > 30) + 1)
#levels(BMI.cat) <- c('Normal', 'Overweight', 'Obese')
#table(BMI.cat)
#boxplot(BMI ~ BMI.cat)
#data.obj0$meta.dat$BMI.cat <- BMI.cat

data.obj0$otu.tab <- data.obj0$abund.list[['Species']]

meta.dat0 <- data.obj0$meta.dat

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw(base_size=16))

study <- 'Castellarin'

dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study))
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/NormalBMI'))
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/CancerBMI'))
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledBMI'))
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS'))

fake.tree <- rcoal(nrow(data.obj0$abund.list[['Species']]))
fake.tree$tip.label <- rownames(data.obj0$abund.list[['Species']])
fake.otu.name <- matrix('Bacteria', nrow(data.obj0$abund.list[['Species']]), 5)
rownames(fake.otu.name) <- rownames(data.obj0$abund.list[['Species']])
phylo.obj0 <- phyloseq(otu_table(data.obj0$abund.list[['Species']], taxa_are_rows=T), phy_tree(fake.tree), 
		tax_table(fake.otu.name), sample_data(data.obj0$meta.dat))



########################################################
# Association with disease status, adjusting for BMI
########################################################
ind <- meta.dat0$Status %in% c('CRC', 'Normal')
data.obj <- subset_data(data.obj0, ind)
dist.obj <- subset_dist(dist.obj0, ind)
data.obj$meta.dat$Status <- factor(data.obj$meta.dat$Status)

phylo.obj <- subset_samples(phylo.obj0, ind)

#### Alpha-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Alpha_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Alpha_Diversity'))


alpha.obj <- perform_alpha_test(data.obj, phylo.obj, rarefy=TRUE, measures=c('Observed',  'Shannon'),  model='lme', 
		random = as.formula(paste0(' ~ 1 | ', 'PatientID')),
		grp.name='Status', ann=paste0(study, '.DS'))

save(alpha.obj,  file=paste0(study, '.alpha.test.RData'))

# Boxplot
generate_alpha_boxplot(data.obj, phylo.obj, rarefy=FALSE, grp.name='Status', 
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2,  7)])', ann=paste0(study, '.Status'))

generate_alpha_boxplot(data.obj, phylo.obj, rarefy=FALSE, grp.name='Status', subject='PatientID',
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2,  7)])', ann=paste0(study, '.StatusPaired'))

#### Beta-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Beta_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Beta_Diversity'))

generate_ordination(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		grp.name='Status',   pca.method='cmd', ann=paste0(study, '.DS'), 
		clab=1.0, cex.pt=1, ellipse=T, cstar= 0, wid=5, hei=5, pdf=TRUE)

generate_ordination(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		grp.name='PatientID',   pca.method='cmd', ann=paste0(study, '.ByPatientID'), 
		clab=0, cex.pt=1, ellipse=T, cstar= 1, wid=5, hei=5, pdf=TRUE)

colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(3, 7)][(x-1) %% 2 + 1]
}
colFnsF <- c(colFunc)
generate_clustering(data.obj, dist.obj, dist.names=c('BC', 'JS'), meta.info='Status', is.labRow=FALSE, ann=paste0(study, '.DS'),  colFnsF=colFnsF)

permanova.obj <- perform_permanova_test(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		PermanovaG.dist=c('BC', 'JS'), strata='PatientID',
		grp.name='Status',  ann=paste0(study, '.DS'))

save(permanova.obj, file=paste0(study, '.beta.test.RData'))
#### Taxa-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Taxa_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Taxa_Diversity'))

# Random slope may be needed
data.obj$meta.dat$PatientID <- factor(data.obj$meta.dat$PatientID )
diff.obj <- perform_differential_analysis_para(data.obj,  grp.name='Status',   RE=TRUE, method='NB', subject='PatientID',
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='TSS', norm.level='Species', intersect.no=4,
		prev=0.1, minp=0.002, medianp=NULL, mt.method='raw', cutoff=0.05, ann=paste0(study, '.DS.TSSNB'))

#visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para', taxa.levels=c('Species'),  ann=paste0(study, '.BMI.OTU.TSSNB'))
visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='Status', test='Para',  subject='PatientID',
		ann=paste0(study, '.DS.taxa.TSSNB'),  colFnsF=colFnsF)
save(diff.obj, file=paste0(study, '.taxa.test.RData'))

