#Raw dada2 output cleanup

#remove eukaryotes, chloroplasts and mitochondria
setwd("C:/Users/jacqu/Desktop/dada2_output")
count_tab <- read.csv("BIOSSCOPE_SizeFrac_seqtab-nochimtaxa_NAchange.csv", header=TRUE,sep=",",row.names=1)
tax_tab <- read.csv("BIOSSCOPE_SizeFrac_taxa_NAchange.csv", header=TRUE,sep=",",row.names=1)
meta <- read.csv("BIOSSCOPE_SizeFrac_metadata.csv", header=TRUE,sep=",",row.names=1)

#Remove rows with chloroplast as the Order name
ASV_nochl <- count_tab %>% filter(Order !="Chloroplast")
TAX_nochl <- tax_tab %>% filter(Order !="Chloroplast")
#Remove rows with mitochondria as the Family name
ASV_nochl.mit <- ASV_nochl %>% filter(Family !="Mitochondria")
TAX_nochl.mit <- TAX_nochl %>% filter(Family !="Mitochondria")
#Remove rows with Eukaryotes as the Kingdom name
ASV_nochl.mit.euk <- ASV_nochl.mit %>% filter(Kingdom !="Eukaryota")
TAX_nochl.mit.euk <- as.matrix(TAX_nochl.mit %>% filter(Kingdom !="Eukaryota"))

#write new files
write.csv(ASV_nochl.mit.euk, file="BIOSSCOPE_SizeFrac_seqtab-nochimtaxa_NoChlMitEuk.csv", row.names=T)
write.csv(TAX_nochl.mit.euk, file="BIOSSCOPE_SizeFrac_taxa_NoChlMitEuk.csv", row.names=T)

#generate rarefaction curve
ASV <- ASV_nochl.mit.euk[,1:(ncol(ASV_nochl.mit.euk)-7)]
rarecurve(ASV, step = 1000, col = "blue", label=F, xlim=c(0,20000))

#rarefy samples to 8000
OTU <- otu_table(ASV, taxa_are_rows = TRUE)
TAX <- tax_table(TAX_nochl.mit.euk)
phy <- phyloseq(OTU,TAX)

set.seed(8800)
rar <- rarefy_even_depth(phy, sample.size = 8000)
ASV_rar <- rar@otu_table
TAX_rar <- rar@tax_table

#calculate relative abundance
ASV_rar_prop <- apply(ASV_rar, 2, function(x) x/sum(x)*100)

#remove singletons (1/8000*100=0.0125%)
single <- rowSums(ASV_rar_prop[,1:ncol(ASV_rar_prop)]) > 0.0125
ASV_NoSingle <- as.data.frame(ASV_rar_prop[single,])

#Remove singletons from rarefied taxonomy file
ASVrow <- as.data.frame(rownames(ASV_NoSingle))
colnames(ASVrow) <- "ASV"
taxa_rar <- as.data.frame(TAX_rar)
taxa_rar$ASV <- rownames(taxa_rar)
taxa_NoSingle <- inner_join(taxa_rarA,ASVrow, by="ASV")

#reverse NA -> NOCHANGE cell contents
taxa_NoSingle$Order <- revalue(taxa_NoSingle$Order, c("NOCHANGE"=NA))
taxa_NoSingle$Family <- revalue(taxa_NoSingle$Family, c("NOCHANGE"=NA))

#create trimmed/rarefied metadata table
sample_rar <- as.data.frame(colnames(ASV_NoSingle))
colnames(sample_rar) <- "Sample"
sample_rar$Sample <- gsub('_','.', sample_rar$Sample)

meta$Sample <- rownames(meta)
meta$Sample <- gsub('_','.', meta$Sample)
meta_all <- full_join(sample_rar,meta, by= "Sample")
meta_rar <- inner_join(sample_rar, meta, by="Sample")

#write new ASV and taxa file
write.csv(taxa_NoSingle, file="BIOSSCOPE_SizeFrac_taxa_rar8000_NoSingle.csv")
write.csv(ASV_NoSingle, file="BIOSSCOPE_SizeFrac_seqtab-nochimtaxa_rar8000_NoSingle.csv")
write.csv(meta_rar,file="BIOSSCOPE_SizeFrac_meta_rar8000.csv")
