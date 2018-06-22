source('~/Dropbox/Workspace/MayoClinic/Stats.R')
load_package()

# 16S data
temp <- load("~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/16S.wk.RData")

process.func <- function (data.obj, tab) {
	cat('Dimension of the functional data: ', dim(tab), '\n')
	cat('Dimension of the taxonomic data: ', dim(data.obj$meta.dat), '\n')
	cat('Number of inconsistent names: ', sum(!(rownames(data.obj$meta.dat) %in% colnames(tab))), '\n')
	
	tab <- tab[, rownames(data.obj$meta.dat)]
	
	cat('Sum of the functional categoris: ', range(colSums(tab)), '\n')
	hist(colSums(tab))
	
	func.names <- rownames(tab)
	func.no <- sapply(strsplit(func.names, ';'), length)
	func.names  <- paste0('L', func.no, '_', sapply(strsplit(func.names, ';'), function (x) x[length(x)]))
	
	cat('No Duplicates:', length(unique(func.names)) == length((func.names)), '\n')
	rownames(tab) <- func.names
	
	abund.list <- data.obj$abund.list
	abund.list[['KEGG_L1']] <- tab[func.no == 1, ]
	abund.list[['KEGG_L2']] <- tab[func.no == 2, ]
	abund.list[['KEGG_L3']] <- tab[func.no == 3, ]
	data.obj$abund.list <- abund.list
	return(data.obj)
}

data.obj <- baxter.stool.data.obj.wk
tab <- as.matrix(read.csv('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/FunctionalData/Baxter_comb1.csv', head = TRUE, row.names = 1))
baxter.stool.data.obj.wk <- process.func(data.obj, tab)

data.obj <- zackular.stool.data.obj.wk
tab <- as.matrix(read.csv('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/FunctionalData/Zackular_comb1.csv', head = TRUE, row.names = 1))
zackular.stool.data.obj.wk <- process.func(data.obj, tab)

data.obj <- zeller.stool.data.obj.wk
tab <- as.matrix(read.csv('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/FunctionalData/Zeller_stool_comb1.csv', head = TRUE, row.names = 1))
zeller.stool.data.obj.wk <- process.func(data.obj, tab)

data.obj <- zeller.tissue.data.obj.wk
tab <- as.matrix(read.csv('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/FunctionalData/Zeller_tissue_comb1.csv', head = TRUE, row.names = 1))
zeller.tissue.data.obj.wk <- process.func(data.obj, tab)


save(baxter.stool.data.obj.wk, baxter.stool.dist.obj.wk, 
     zackular.stool.data.obj.wk, zackular.stool.dist.obj.wk,
     zeller.stool.data.obj.wk, zeller.stool.dist.obj.wk,  
     zeller.tissue.data.obj.wk, zeller.tissue.dist.obj.wk, file = '~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/16S.wk2.RData')

# WGS data
temp <- load("~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/WGS.wk.RData")

process.func <- function (data.obj, tab) {

	cat('Dimension of the functional data: ', dim(tab), '\n')
	cat('Dimension of the taxonomic data: ', dim(data.obj$meta.dat), '\n')
	cat('Number of inconsistent names: ', sum(!(rownames(data.obj$meta.dat) %in% colnames(tab))), '\n')
	
	samIDs <- intersect(rownames(data.obj$meta.dat), colnames(tab))
	data.obj <- subset_data(data.obj, samIDs)
	
	tab <- tab[, samIDs]
	
	cat('Sum of the functional categoris: ', range(colSums(tab)), '\n')
	hist(colSums(tab))
	
	func.names <- rownames(tab)
    func.ind <- grepl('^k', func.names)
	
	abund.list <- data.obj$abund.list
	abund.list[['KEGG']] <- tab[func.ind, ]
	abund.list[['Module']] <- tab[!func.ind == 2, ]
	
	data.obj$abund.list <- abund.list
	return(data.obj)
}

data.obj <- Castellarin.data.obj.wk
tab <- as.matrix(read.csv('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/FunctionalData/Castellarin_comb1.csv', head = TRUE, row.names = 1))
Castellarin.data.obj.wk <- process.func(data.obj, tab)
Castellarin..dist.obj.wk <- subset_dist(Castellarin.dist.obj.wk, rownames(Castellarin.data.obj.wk$meta.dat))

# Three inconsistent names - we lose three samples
data.obj <- Feng.data.obj.wk
tab <- as.matrix(read.csv('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/FunctionalData/Feng_comb1.csv', head = TRUE, row.names = 1))
Feng.data.obj.wk <- process.func(data.obj, tab)
Feng.dist.obj.wk <- subset_dist(Feng.dist.obj.wk, rownames(Feng.data.obj.wk$meta.dat))

data.obj <- Voghtmann.data.obj.wk
tab <- as.matrix(read.csv('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/FunctionalData/Voghtmann_comb1.csv', head = TRUE, row.names = 1))
Voghtmann.data.obj.wk <- process.func(data.obj, tab)
Voghtmann.dist.obj.wk <- subset_dist(Voghtmann.dist.obj.wk, rownames(Voghtmann.data.obj.wk$meta.dat))

data.obj <- Zeller.data.obj.wk
tab <- as.matrix(read.csv('~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/FunctionalData/Zeller_comb1.csv', head = TRUE, row.names = 1))
Zeller.data.obj.wk <- process.func(data.obj, tab)
Zeller.dist.obj.wk <- subset_dist(Zeller.dist.obj.wk, rownames(Zeller.data.obj.wk$meta.dat))

save(Castellarin.data.obj.wk, Castellarin.dist.obj.wk,
     Feng.data.obj.wk,        Feng.dist.obj.wk,       
     Voghtmann.data.obj.wk,   Voghtmann.dist.obj.wk,  
     Zeller.data.obj.wk,      Zeller.dist.obj.wk, file = '~/Dropbox/Workspace/MayoClinic/Collaboration/2016_11_16_Leigh_Metaanalysis/Data/WGS.wk2.RData')

