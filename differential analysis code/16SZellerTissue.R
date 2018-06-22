source('~/Dropbox/Workspace/MayoClinic/Stats.R')
load_package()

load("~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Data/16S.wk.RData")

data.obj0 <- zeller.tissue.data.obj.wk
dist.obj0 <- zeller.tissue.dist.obj.wk

BMI <- data.obj0$meta.dat$BMI
BMI.cat <- factor((BMI > 25) + (BMI > 30) + 1)
levels(BMI.cat) <- c('Normal', 'Overweight', 'Obese')
table(BMI.cat)
boxplot(BMI ~ BMI.cat)
data.obj0$meta.dat$BMI.cat <- BMI.cat

meta.dat0 <- data.obj0$meta.dat

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw(base_size=16))

study <- 'zeller.tissue'

dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study))
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/NormalBMI'))
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/CancerBMI'))
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/PooledDS'))

########################################################
# Association in normal patients only
########################################################
ind <- meta.dat0$Status == 'Normal'
data.obj <- subset_data(data.obj0, ind)
dist.obj <- subset_dist(dist.obj0, ind)
phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
		tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))

#### Alpha-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/NormalBMI/Alpha_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/NormalBMI/Alpha_Diversity'))

min(colSums(data.obj$otu.tab))
alpha.obj.cat <- perform_alpha_test(data.obj, phylo.obj, rarefy=FALSE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='BMI.cat', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI.cat'))

alpha.obj <- perform_alpha_test(data.obj, phylo.obj, rarefy=FALSE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='BMI', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI'))

save(alpha.obj, alpha.obj.cat, file=paste0(study, '.alpha.test.RData'))

# Scatter plot
BMI <- data.obj$meta.dat$BMI
alpha.diversity <- alpha.obj$alpha.diversity

colnames(alpha.diversity)
df <- data.frame(BMI=BMI, alpha.diversity)
df <- melt(df, id.vars=c('BMI'))
df$variable <- factor(df$variable)
levels(df$variable) <- c('Observed OTU number', 'Shannon index')
ggplot(df, aes(BMI, value)) +
		geom_point(col="#0072B2", alpha=0.5, size=2) +
		geom_smooth(method='lm') +
		facet_wrap(~variable, ncol=2, scales='free_y') +
		theme_bw()
ggsave(paste0('Alpha_diversity_Scatterplot_', study, '.BMI.pdf'), width=7, height=4)

# Boxplot
generate_alpha_boxplot(data.obj, phylo.obj, rarefy=FALSE, grp.name='BMI.cat', 
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2, 4, 7)])', ann=paste0(study, '.BMI.cat'))

#### Beta-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/NormalBMI/Beta_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/NormalBMI/Beta_Diversity'))

generate_ordination(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name='BMI.cat',   pca.method='cmd', ann=paste0(study, '.BMI.cat'), 
		clab=1.0, cex.pt=1, ellipse=T, cstar= 0, wid=5, hei=5, pdf=TRUE)

colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(2, 4, 7)][(x-1) %% 3 + 1]
}
colFnsF <- c(colFunc)
generate_clustering(data.obj, dist.obj, meta.info='BMI.cat', is.labRow=FALSE, ann=paste0(study, '.BMI.cat'),  colFnsF=colFnsF)

permanova.obj <- perform_permanova_test(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		PermanovaG.dist=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name='BMI', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI'))
save(permanova.obj, file=paste0(study, '.beta.test.RData'))

#### Taxa-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/NormalBMI/Taxa_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/NormalBMI/Taxa_Diversity'))

diff.obj <- perform_differential_analysis_para(data.obj,  grp.name='BMI', adj.name=c('Gender', 'Age'),  RE=FALSE, method='NB', 
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='TSS', norm.level='Species', intersect.no=4,
		prev=0.1, minp=0.002, medianp=NULL, mt.method='raw', cutoff=0.05, ann=paste0(study, '.BMI.TSSNB'))

#visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para', taxa.levels=c('Species'),  ann=paste0(study, '.BMI.OTU.TSSNB'))
visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para',  ann=paste0(study, '.BMI.taxa.TSSNB'),  colFnsF=colFnsF)

save(diff.obj, file=paste0(study, '.taxa.test.RData'))

########################################################
# Association in cancer patients only
########################################################
ind <- meta.dat0$Status == 'CRC'
data.obj <- subset_data(data.obj0, ind)
dist.obj <- subset_dist(dist.obj0, ind)
phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
		tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))

#### Alpha-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/CancerBMI/Alpha_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/CancerBMI/Alpha_Diversity'))

min(colSums(data.obj$otu.tab))
alpha.obj.cat <- perform_alpha_test(data.obj, phylo.obj, rarefy=FALSE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='BMI.cat', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI.cat'))

alpha.obj <- perform_alpha_test(data.obj, phylo.obj, rarefy=FALSE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='BMI', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI'))

save(alpha.obj, alpha.obj.cat, file=paste0(study, '.alpha.test.RData'))

# Scatter plot
BMI <- data.obj$meta.dat$BMI
alpha.diversity <- alpha.obj$alpha.diversity

colnames(alpha.diversity)
df <- data.frame(BMI=BMI, alpha.diversity)
df <- melt(df, id.vars=c('BMI'))
df$variable <- factor(df$variable)
levels(df$variable) <- c('Observed OTU number', 'Shannon index')
ggplot(df, aes(BMI, value)) +
		geom_point(col="#0072B2", alpha=0.5, size=2) +
		geom_smooth(method='lm') +
		facet_wrap(~variable, ncol=2, scales='free_y') +
		theme_bw()
ggsave(paste0('Alpha_diversity_Scatterplot_', study, '.BMI.pdf'), width=7, height=4)

# Boxplot
generate_alpha_boxplot(data.obj, phylo.obj, rarefy=FALSE, grp.name='BMI.cat', 
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2, 4, 7)])', ann=paste0(study, '.BMI.cat'))

#### Beta-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/CancerBMI/Beta_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/CancerBMI/Beta_Diversity'))

generate_ordination(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name='BMI.cat',   pca.method='cmd', ann=paste0(study, '.BMI.cat'), 
		clab=1.0, cex.pt=1, ellipse=T, cstar= 0, wid=5, hei=5, pdf=TRUE)

colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(2, 4, 7)][(x-1) %% 3 + 1]
}
colFnsF <- c(colFunc)
generate_clustering(data.obj, dist.obj, meta.info='BMI.cat', is.labRow=FALSE, ann=paste0(study, '.BMI.cat'),  colFnsF=colFnsF)

permanova.obj <- perform_permanova_test(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		PermanovaG.dist=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name='BMI', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI'))

save(permanova.obj, file=paste0(study, '.beta.test.RData'))
#### Taxa-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/CancerBMI/Taxa_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/CancerBMI/Taxa_Diversity'))

diff.obj <- perform_differential_analysis_para(data.obj,  grp.name='BMI', adj.name=c('Gender', 'Age'),  RE=FALSE, method='NB', 
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='TSS', norm.level='Species', intersect.no=4,
		prev=0.1, minp=0.002, medianp=NULL, mt.method='raw', cutoff=0.05, ann=paste0(study, '.BMI.TSSNB'))

#visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para', taxa.levels=c('Species'),  ann=paste0(study, '.BMI.OTU.TSSNB'))
visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para',  ann=paste0(study, '.BMI.taxa.TSSNB'),  colFnsF=colFnsF)
save(diff.obj, file=paste0(study, '.taxa.test.RData'))


########################################################
# Association with disease status, adjusting for BMI
########################################################
ind <- meta.dat0$Status %in% c('CRC', 'Normal')
data.obj <- subset_data(data.obj0, ind)
dist.obj <- subset_dist(dist.obj0, ind)
data.obj$meta.dat$Status <- factor(data.obj$meta.dat$Status)

phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
		tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))

#### Alpha-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/PooledDS/Alpha_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/PooledDS/Alpha_Diversity'))

min(colSums(data.obj$otu.tab))

alpha.obj <- perform_alpha_test(data.obj, phylo.obj, rarefy=FALSE, measures=c('Observed',  'Shannon'),  model='lme', 
		random = as.formula(paste0(' ~ 1 | ', 'PatientID')),
		grp.name='Status', adj.name=c('Gender', 'Age', 'BMI'), ann=paste0(study, '.DS'))

save(alpha.obj,  file=paste0(study, '.alpha.test.RData'))

# Boxplot
generate_alpha_boxplot(data.obj, phylo.obj, rarefy=FALSE, grp.name='Status', 
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2,  7)])', ann=paste0(study, '.Status'))

generate_alpha_boxplot(data.obj, phylo.obj, rarefy=FALSE, grp.name='Status', subject='PatientID',
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2,  7)])', ann=paste0(study, '.StatusPaired'))

#### Beta-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/PooledDS/Beta_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/PooledDS/Beta_Diversity'))

generate_ordination(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name='Status',   pca.method='cmd', ann=paste0(study, '.DS'), 
		clab=1.0, cex.pt=1, ellipse=T, cstar= 0, wid=5, hei=5, pdf=TRUE)

generate_ordination(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name='PatientID',   pca.method='cmd', ann=paste0(study, '.ByPatientID'), 
		clab=0, cex.pt=1, ellipse=T, cstar= 1, wid=5, hei=5, pdf=TRUE)

colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(3, 7)][(x-1) %% 2 + 1]
}
colFnsF <- c(colFunc)
generate_clustering(data.obj, dist.obj, meta.info='Status', is.labRow=FALSE, ann=paste0(study, '.DS'),  colFnsF=colFnsF)

permanova.obj <- perform_permanova_test(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		PermanovaG.dist=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), strata='PatientID',
		grp.name='Status', adj.name=c('Gender', 'Age', 'BMI'), ann=paste0(study, '.DS'))

save(permanova.obj, file=paste0(study, '.beta.test.RData'))
#### Taxa-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/PooledDS/Taxa_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/16S/', study, '/PooledDS/Taxa_Diversity'))

data.obj$meta.dat$PatientID <- factor(data.obj$meta.dat$PatientID )

diff.obj <- perform_differential_analysis_para(data.obj,  grp.name='Status', adj.name=c('Gender', 'Age', 'BMI'),  RE=TRUE, method='NB', subject='PatientID',
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='TSS', norm.level='Species', intersect.no=4,
		prev=0.1, minp=0.002, medianp=NULL, mt.method='raw', cutoff=0.05, ann=paste0(study, '.DS.TSSNB'))

#visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para', taxa.levels=c('Species'),  ann=paste0(study, '.BMI.OTU.TSSNB'))
visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='Status', test='Para',  ann=paste0(study, '.DS.taxa.TSSNB'),  subject='PatientID', colFnsF=colFnsF)
visualize_differential_analysis(data.obj, diff.obj,  mt.method='fdr', cutoff=0.05, grp.name='Status', test='Para',  ann=paste0(study, '.DS.taxa.TSSNB'),  subject='PatientID', colFnsF=colFnsF)
save(diff.obj, file=paste0(study, '.taxa.test.RData'))

