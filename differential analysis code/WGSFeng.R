source('~/Dropbox/Workspace/MayoClinic/Stats.R')
load_package()

load("~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Data/WGS.wk.RData")

data.obj0 <- Feng.data.obj.wk
dist.obj0 <- Feng.dist.obj.wk

BMI <- data.obj0$meta.dat$BMI
BMI.cat <- factor((BMI > 25) + (BMI > 30) + 1)
levels(BMI.cat) <- c('Normal', 'Overweight', 'Obese')
table(BMI.cat)
boxplot(BMI ~ BMI.cat)
data.obj0$meta.dat$BMI.cat <- BMI.cat

data.obj0$otu.tab <- data.obj0$abund.list[['Species']]

meta.dat0 <- data.obj0$meta.dat

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw(base_size=16))

study <- 'Feng'

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
# Association in normal patients only
########################################################
ind <- meta.dat0$Status == 'Normal'
data.obj <- subset_data(data.obj0, ind)
dist.obj <- subset_dist(dist.obj0, ind)

phylo.obj <- subset_samples(phylo.obj0, ind)

#### Alpha-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/NormalBMI/Alpha_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/NormalBMI/Alpha_Diversity'))

alpha.obj.cat <- perform_alpha_test(data.obj, phylo.obj, rarefy=TRUE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='BMI.cat', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI.cat'))

alpha.obj <- perform_alpha_test(data.obj, phylo.obj, rarefy=TRUE, measures=c('Observed',  'Shannon'),  model='lm', 
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
generate_alpha_boxplot(data.obj, phylo.obj, rarefy=TRUE, grp.name='BMI.cat', 
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2, 4, 7)])', ann=paste0(study, '.BMI.cat'))

#### Beta-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/NormalBMI/Beta_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/NormalBMI/Beta_Diversity'))

generate_ordination(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		grp.name='BMI.cat',   pca.method='cmd', ann=paste0(study, '.BMI.cat'), 
		clab=1.0, cex.pt=1, ellipse=T, cstar= 0, wid=5, hei=5, pdf=TRUE)

colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(2, 4, 7)][(x-1) %% 3 + 1]
}
colFnsF <- c(colFunc)
generate_clustering(data.obj, dist.obj, dist.names=c('BC', 'JS'), meta.info='BMI.cat', is.labRow=FALSE, ann=paste0(study, '.BMI.cat'),  colFnsF=colFnsF)

permanova.obj <- perform_permanova_test(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		PermanovaG.dist=c('BC', 'JS'),
		grp.name='BMI', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI'))
save(permanova.obj, file=paste0(study, '.beta.test.RData'))

#### Taxa-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/NormalBMI/Taxa_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/NormalBMI/Taxa_Diversity'))

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

phylo.obj <- subset_samples(phylo.obj0, ind)

#### Alpha-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/CancerBMI/Alpha_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/CancerBMI/Alpha_Diversity'))


alpha.obj.cat <- perform_alpha_test(data.obj, phylo.obj, rarefy=TRUE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='BMI.cat', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI.cat'))

alpha.obj <- perform_alpha_test(data.obj, phylo.obj, rarefy=TRUE, measures=c('Observed',  'Shannon'),  model='lm', 
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
generate_alpha_boxplot(data.obj, phylo.obj, rarefy=TRUE, grp.name='BMI.cat', 
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2, 4, 7)])', ann=paste0(study, '.BMI.cat'))

#### Beta-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/CancerBMI/Beta_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/CancerBMI/Beta_Diversity'))

generate_ordination(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		grp.name='BMI.cat',   pca.method='cmd', ann=paste0(study, '.BMI.cat'), 
		clab=1.0, cex.pt=1, ellipse=T, cstar= 0, wid=5, hei=5, pdf=TRUE)

colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(2, 4, 7)][(x-1) %% 3 + 1]
}
colFnsF <- c(colFunc)
generate_clustering(data.obj, dist.obj, dist.names=c('BC', 'JS'), meta.info='BMI.cat', is.labRow=FALSE, ann=paste0(study, '.BMI.cat'),  colFnsF=colFnsF)

permanova.obj <- perform_permanova_test(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		PermanovaG.dist=c('BC', 'JS'),
		grp.name='BMI', adj.name=c('Gender', 'Age'), ann=paste0(study, '.BMI'))

save(permanova.obj, file=paste0(study, '.beta.test.RData'))
#### Taxa-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/CancerBMI/Taxa_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/CancerBMI/Taxa_Diversity'))

diff.obj <- perform_differential_analysis_para(data.obj,  grp.name='BMI', adj.name=c('Gender', 'Age'),  RE=FALSE, method='NB', 
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='TSS', norm.level='Species', intersect.no=4,
		prev=0.1, minp=0.002, medianp=NULL, mt.method='raw', cutoff=0.05, ann=paste0(study, '.BMI.TSSNB'))

#visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para', taxa.levels=c('Species'),  ann=paste0(study, '.BMI.OTU.TSSNB'))
visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para',  ann=paste0(study, '.BMI.taxa.TSSNB'),  colFnsF=colFnsF)
save(diff.obj, file=paste0(study, '.taxa.test.RData'))

########################################################
# Association in pooled patients only, adjusting for disease status
########################################################
ind <- meta.dat0$Status %in% c('CRC', 'Normal', 'Adenoma')
data.obj <- subset_data(data.obj0, ind)
dist.obj <- subset_dist(dist.obj0, ind)

phylo.obj <- subset_samples(phylo.obj0, ind)

#### Alpha-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledBMI/Alpha_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledBMI/Alpha_Diversity'))


alpha.obj.cat <- perform_alpha_test(data.obj, phylo.obj, rarefy=TRUE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='BMI.cat', adj.name=c('Gender', 'Age', 'Status'), ann=paste0(study, '.BMI.cat'))

alpha.obj <- perform_alpha_test(data.obj, phylo.obj, rarefy=TRUE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='BMI', adj.name=c('Gender', 'Age', 'Status'), ann=paste0(study, '.BMI'))

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
generate_alpha_boxplot(data.obj, phylo.obj, rarefy=TRUE, grp.name='BMI.cat', 
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2, 4, 7)])', ann=paste0(study, '.BMI.cat'))

#### Beta-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledBMI/Beta_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledBMI/Beta_Diversity'))

generate_ordination(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		grp.name='BMI.cat',   pca.method='cmd', ann=paste0(study, '.BMI.cat'), 
		clab=1.0, cex.pt=1, ellipse=T, cstar= 0, wid=5, hei=5, pdf=TRUE)

colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(2, 4, 7)][(x-1) %% 3 + 1]
}
colFnsF <- c(colFunc)
generate_clustering(data.obj, dist.obj, dist.names=c('BC', 'JS'), meta.info='BMI.cat', is.labRow=FALSE, ann=paste0(study, '.BMI.cat'),  colFnsF=colFnsF)

permanova.obj <- perform_permanova_test(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		PermanovaG.dist=c('BC', 'JS'),
		grp.name='BMI', adj.name=c('Gender', 'Age', 'Status'), ann=paste0(study, '.BMI'))

save(permanova.obj, file=paste0(study, '.beta.test.RData'))
#### Taxa-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledBMI/Taxa_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledBMI/Taxa_Diversity'))

diff.obj <- perform_differential_analysis_para(data.obj,  grp.name='BMI', adj.name=c('Gender', 'Age', 'Status'),  RE=FALSE, method='NB', 
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

phylo.obj <- subset_samples(phylo.obj0, ind)

#### Alpha-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Alpha_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Alpha_Diversity'))


alpha.obj <- perform_alpha_test(data.obj, phylo.obj, rarefy=TRUE, measures=c('Observed',  'Shannon'),  model='lm', 
		grp.name='Status', adj.name=c('Gender', 'Age', 'BMI'), ann=paste0(study, '.DS'))

save(alpha.obj,  file=paste0(study, '.alpha.test.RData'))

# Boxplot
generate_alpha_boxplot(data.obj, phylo.obj, rarefy=FALSE, grp.name='Status', 
		measures=c('Observed', 'Shannon'), gg.cmd='scale_colour_manual(values=cbPalette[c(2,  7)])', ann=paste0(study, '.Status'))

#### Beta-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Beta_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Beta_Diversity'))

generate_ordination(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		grp.name='Status',   pca.method='cmd', ann=paste0(study, '.DS'), 
		clab=1.0, cex.pt=1, ellipse=T, cstar= 0, wid=5, hei=5, pdf=TRUE)

colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(3, 7)][(x-1) %% 2 + 1]
}
colFnsF <- c(colFunc)
generate_clustering(data.obj, dist.obj, dist.names=c('BC', 'JS'), meta.info='Status', is.labRow=FALSE, ann=paste0(study, '.DS'),  colFnsF=colFnsF)

permanova.obj <- perform_permanova_test(data.obj, dist.obj, dist.names=c('BC', 'JS'), 
		PermanovaG.dist=c('BC', 'JS'),
		grp.name='Status', adj.name=c('Gender', 'Age', 'BMI'), ann=paste0(study, '.DS'))

save(permanova.obj, file=paste0(study, '.beta.test.RData'))
#### Taxa-diversity comparison
dir.create(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Taxa_Diversity'))
setwd(paste0('~/Dropbox/Workspace/MayoClinic/2016_11_16_Leigh_Metaanalysis/Result/WGS/', study, '/PooledDS/Taxa_Diversity'))

diff.obj <- perform_differential_analysis_para(data.obj,  grp.name='Status', adj.name=c('Gender', 'Age', 'BMI'),  RE=FALSE, method='NB', 
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='TSS', norm.level='Species', intersect.no=4,
		prev=0.1, minp=0.002, medianp=NULL, mt.method='raw', cutoff=0.05, ann=paste0(study, '.DS.TSSNB'))

#visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat', test='Para', taxa.levels=c('Species'),  ann=paste0(study, '.BMI.OTU.TSSNB'))
visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='Status', test='Para',  ann=paste0(study, '.DS.taxa.TSSNB'),  colFnsF=colFnsF)
save(diff.obj, file=paste0(study, '.taxa.test.RData'))

