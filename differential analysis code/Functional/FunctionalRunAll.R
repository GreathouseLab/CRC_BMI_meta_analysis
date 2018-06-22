

source('~/Dropbox/Workspace/MayoClinic/Stats.R')
load_package()

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colFunc <- function (x) {
	x <- 1:x
	cbPalette[c(2, 4, 7)][(x-1) %% 3 + 1]
}
colFnsF <- c(colFunc)
theme_set(theme_bw(base_size=16))
###########################################################################################
cat('Run 16S ...\n')
load("~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/16S.wk2.RData")

for (study in c('baxter.stool', 'zackular.stool', 'zeller.stool', 'zeller.tissue')) {
	
	cat(study, '...\n')
	data.obj0 <- get(paste0(study, ".data.obj.wk"))
	dist.obj0 <- get(paste0(study, ".dist.obj.wk"))
	
	adj.name <- c('Gender', 'Age')
	
	BMI <- data.obj0$meta.dat$BMI
	BMI.cat <- factor((BMI > 25) + (BMI > 30) + 1)
	levels(BMI.cat) <- c('Normal', 'Overweight', 'Obese')
	table(BMI.cat)
	boxplot(BMI ~ BMI.cat)
	data.obj0$meta.dat$BMI.cat <- BMI.cat
	
	meta.dat0 <- data.obj0$meta.dat
	
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study))
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/NormalBMI'))
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/CancerBMI'))

	########################################################
    # Association in normal patients only
	########################################################
	ind <- meta.dat0$Status == 'Normal'
	data.obj <- subset_data(data.obj0, ind)
	dist.obj <- subset_dist(dist.obj0, ind)
	
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/NormalBMI'))
	setwd(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/NormalBMI'))
	
	diff.obj <- perform_differential_analysis(data.obj,  grp.name='BMI', adj.name=adj.name, method='perm', 
			taxa.levels=c('KEGG_L1', 'KEGG_L2', 'KEGG_L3'), norm='TSS', norm.level='KEGG_L1', 
			prev=0.1, minp=0.001, medianp=NULL, cv = 0.05, mt.method='raw', cutoff=0.05, ann=paste0(study, '.BMI.TSS.PERM'))
	
	visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat',  taxa.levels = c('KEGG_L1', 'KEGG_L2', 'KEGG_L3'), 
			ann=paste0(study, '.BMI.TSS.PERM'),  colFnsF=colFnsF)
	
	save(diff.obj, file=paste0(study, '.func.test.RData'))
	
	########################################################
    # Association in cancer patients only
	########################################################
	ind <- meta.dat0$Status == 'CRC'
	data.obj <- subset_data(data.obj0, ind)
	dist.obj <- subset_dist(dist.obj0, ind)
	
	#### Taxa-diversity comparison
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/CancerBMI'))
	setwd(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/CancerBMI'))
	
	diff.obj <- perform_differential_analysis(data.obj,  grp.name='BMI', adj.name=adj.name, method='perm', 
			taxa.levels=c('KEGG_L1', 'KEGG_L2', 'KEGG_L3'), norm='TSS', norm.level='KEGG_L1', 
			prev=0.1, minp=0.001, medianp=NULL, cv = 0.05, mt.method='raw', cutoff=0.05, ann=paste0(study, '.BMI.TSS.PERM'))

	visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat',  taxa.levels = c('KEGG_L1', 'KEGG_L2', 'KEGG_L3'), 
			ann=paste0(study, '.BMI.TSS.PERM'),  colFnsF=colFnsF)
	
	save(diff.obj, file=paste0(study, '.func.test.RData'))
	
}

###########################################################################################
cat('Run WGS ...\n')
temp <- load("~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/WGS.wk2.RData")

for (study in c('Voghtmann', 'Feng', 'Zeller')) {
	
	cat(study, '...\n')
	data.obj0 <- get(paste0(study, ".data.obj.wk"))
	dist.obj0 <- get(paste0(study, ".dist.obj.wk"))
	
	if (study == 'Zeller') {
		adj.name <- c('Gender', 'Age', 'Country')
	} else {
		adj.name <- c('Gender', 'Age')
	}
	
	
	BMI <- data.obj0$meta.dat$BMI
	BMI.cat <- factor((BMI > 25) + (BMI > 30) + 1)
	levels(BMI.cat) <- c('Normal', 'Overweight', 'Obese')
	table(BMI.cat)
	boxplot(BMI ~ BMI.cat)
	data.obj0$meta.dat$BMI.cat <- BMI.cat
	data.obj0$otu.tab <- data.obj0$abund.list[['Species']] # Place holder
	
	meta.dat0 <- data.obj0$meta.dat
	
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study))
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/NormalBMI'))
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/CancerBMI'))
	
	########################################################
	# Association in normal patients only
	########################################################
	ind <- meta.dat0$Status == 'Normal'
	data.obj <- subset_data(data.obj0, ind)
	dist.obj <- subset_dist(dist.obj0, ind)
	
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/NormalBMI'))
	setwd(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/NormalBMI'))
	
	diff.obj <- perform_differential_analysis(data.obj,  grp.name='BMI', adj.name=adj.name, method='perm', 
			taxa.levels=c('KEGG', 'Module'), norm='TSS', norm.level='KEGG', 
			prev=0.1, minp=0.001, medianp=NULL, cv = 0.05, mt.method='raw', cutoff=0.05, ann=paste0(study, '.BMI.TSS.PERM'))
	
	visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat',  taxa.levels = c('KEGG', 'Module'), 
			ann=paste0(study, '.BMI.TSS.PERM'),  colFnsF=colFnsF)
	
	save(diff.obj, file=paste0(study, '.func.test.RData'))
	
	########################################################
	# Association in cancer patients only
	########################################################
	ind <- meta.dat0$Status == 'CRC'
	data.obj <- subset_data(data.obj0, ind)
	dist.obj <- subset_dist(dist.obj0, ind)
	
	#### Taxa-diversity comparison
	dir.create(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/CancerBMI'))
	setwd(paste0('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Result/Function/', study, '/CancerBMI'))
	
	diff.obj <- perform_differential_analysis(data.obj,  grp.name='BMI', adj.name=adj.name, method='perm', 
			taxa.levels=c('KEGG', 'Module'), norm='TSS', norm.level='KEGG', 
			prev=0.1, minp=0.001, medianp=NULL, cv = 0.05, mt.method='raw', cutoff=0.05, ann=paste0(study, '.BMI.TSS.PERM'))
	
	visualize_differential_analysis(data.obj, diff.obj,  mt.method='raw', cutoff=0.05, grp.name='BMI.cat',  taxa.levels = c('KEGG', 'Module'), 
			ann=paste0(study, '.BMI.TSS.PERM'),  colFnsF=colFnsF)
	
	save(diff.obj, file=paste0(study, '.func.test.RData'))
}





