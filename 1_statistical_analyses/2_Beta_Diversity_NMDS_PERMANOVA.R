#Run NMDS ordination analyses and PERMANOVAs
#import data
ASV <- read.csv("BIOSSCOPE_SizeFrac_seqtab-nochimtaxa_rar8000_NoSingle.csv", header=TRUE,sep=",",row.names=1)
taxa <- read.csv("BIOSSCOPE_SizeFrac_taxa_rar8000_NoSingle.csv", header=TRUE,sep=",",row.names=1, na.strings = "NOPE")
meta <- read.csv("BIOSSCOPE_SizeFrac_meta_rar8000.csv", header=TRUE,sep=",",row.names=1)

#rename ASV names in ASV table from DNA string to taxonomy
taxa_num <- as.data.frame(c(1:nrow(taxa)))
taxa_concat <- str_c(taxa$Kingdom, taxa$Phylum, taxa$Class, taxa$Order, taxa$Family, taxa$Genus, taxa$Species,  taxa_num$`c(1:nrow(taxa))`, sep="_")
t_ASV <- as.data.frame(t(ASV))
colnames(t_ASV) <- taxa_concat

#wrangle the different subsets
t_ASV$Sample <- rownames(t_ASV)
t_ASV$Sample <- gsub('_','.', t_ASV$Sample)
mydata <- inner_join(meta, t_ASV, by="Sample")
rownames(mydata) <- mydata$Sample
mydata$Fraction_size_um <- as.character(mydata$Fraction_size_um)
mydata$Depth_m <- as.numeric(mydata$Depth_m)

mydata <- as.data.frame(read.csv("BIOSSCOPE_SizeFrac_EnvASV_noBLANKv3.csv", row.names = 1))
mydata$Fraction_size_um <- as.character(mydata$Fraction_size_um)
all.asin <- as.matrix(asin(sqrt(mydata[,9:ncol(mydata)]/100)))

#ordinate ALL samples
ASV.bc<-as.matrix(vegdist(ASV_asin, method="bray"))
ALLnmds <- metaMDS(all.asin, distance = "bray")
ALL.scores = as.data.frame(scores(ALLnmds))
write.csv(ALL.scores,"BIOSSCOPE_SizeFrac_ALL_NMDSscores.csv")

  #plot with depth as color
ggplot(data = ALL.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.8, aes(shape=mydata$Fraction_size_um, color=mydata$Depth_m))+ 
        ggtitle("NMDS ordination of all pump samples")+
        theme_classic(base_size = 14)+ 
        scale_color_viridis(option = "C", direction = -1)

  #plot with depth as color but with different colors
pal <- rev(wes_palette("Zissou1", 7, type = "continuous"))
ggplot(data = ALL.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.8, aes(shape=mydata$Fraction_size_um, color = mydata$depth.bin.fine))+ 
        ggtitle("NMDS ordination of all pump samples")+
        theme_classic(base_size = 14)+
        scale_color_manual(values = pal)

  #plot with cruise as color
ggplot(data = ALL.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.4, aes(shape=mydata$Fraction_size_um, color=mydata$Cruise))+ 
        ggtitle("NMDS ordination of all pump samples")+
        theme_classic(base_size = 14)

#ordinate samples by depth bin (code to generate asin files are below in alpha diversity section)
UE.nmds <- metaMDS(UE.asin, distance = "bray")
UE.scores <- as.data.frame(scores(UE.nmds))
DCM.nmds <- metaMDS(DCM.asin, distance = "bray")
DCM.scores <- as.data.frame(scores(DCM.nmds))
MESO.nmds <- metaMDS(MESO.asin, distance = "bray")
MESO.scores <- as.data.frame(scores(MESO.nmds))

ggplot(data = UE.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.4, aes(shape=UE$Cruise, color=UE$Fraction_size_um))+ 
        ggtitle("NMDS ordination of UE")+
        theme_classic(base_size = 14)

ggplot(data = DCM.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.4, aes(shape=DCM$Cruise, color=DCM$Fraction_size_um))+ 
        ggtitle("NMDS ordination of DCM")+
        theme_classic(base_size = 14)

ggplot(data = MESO.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.4, aes(shape=MESO$Cruise, color=MESO$Fraction_size_um))+ 
        ggtitle("NMDS ordination of MESO")+
        theme_classic(base_size = 14)

#ordinate by fraction
frac0.2.nmds <- metaMDS(frac0.2.asin, distance = "bray")
frac0.2.scores <- as.data.frame(scores(frac0.2.nmds))
frac1.2.nmds <- metaMDS(frac1.2.asin, distance = "bray")
frac1.2.scores <- as.data.frame(scores(frac1.2.nmds))
frac5.nmds <- metaMDS(frac5.asin, distance = "bray")
frac5.scores <- as.data.frame(scores(frac5.nmds))
frac20.nmds <- metaMDS(frac20.asin, distance = "bray")
frac20.scores <- as.data.frame(scores(frac20.nmds))

ggplot(data = frac0.2.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.4, aes(shape=frac0.2$Cruise, color=frac0.2$depth.bin))+ 
        ggtitle("NMDS ordination of 0.2um fraction")+
        theme_classic(base_size = 14)
ggplot(data = frac1.2.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.4, aes(shape=frac1.2$Cruise, color=frac1.2$depth.bin))+ 
        ggtitle("NMDS ordination of 1.2um fraction")+
        theme_classic(base_size = 14)
ggplot(data = frac5.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.4, aes(shape=frac5$Cruise, color=frac5$depth.bin))+ 
        ggtitle("NMDS ordination of 5um fraction")+
        theme_classic(base_size = 14)
ggplot(data = frac20.scores, 
       mapping = aes(x=NMDS1, y=NMDS2))+ 
        geom_point(size=3, alpha=0.4, aes(shape=frac20$Cruise, color=frac20$depth.bin))+ 
        ggtitle("NMDS ordination of 20um fraction")+
        theme_classic(base_size = 14)

#Run PERMANOVAs
All.adonis1 <- adonis(ASV.bc ~ Cruise, data = mydata)
All.adonis2 <- adonis(ASV.bc ~ Depth_m, data = mydata)
All.adonis3 <- adonis(ASV.bc ~ Fraction_size_um, data = mydata)

UE.adonis <- adonis(UE.asin ~ Fraction_size_um, data = UE, method = "bray")
DCM.adonis <- adonis(DCM.asin ~ Fraction_size_um, data = DCM, method = "bray")
MESO.adonis <- adonis(MESO.asin ~ Fraction_size_um, data = MESO, method = "bray")

frac0.2.adonis <- adonis(frac0.2.asin ~ depth.bin, data = frac0.2, method = "bray")
frac1.2.adonis <- adonis(frac1.2.asin ~ depth.bin, data = frac1.2, method = "bray")
frac5.adonis <- adonis(frac5.asin ~ depth.bin, data = frac5, method = "bray")
frac20.adonis <- adonis(frac20.asin ~ depth.bin, data = frac20, method = "bray")

#Run pairwise PERMANOVAs
All.pair1 <- pairwise.adonis(ASV.bc,mydata$Cruise)
All.pair2 <- pairwise.adonis(ASV.bc,mydata$Depth_m)
All.pair3 <- pairwise.adonis(ASV.bc,mydata$Fraction_size_um)
All.pair4 <- pairwise.adonis(ASV_asin,mydata$fraction.depth.bin, sim.method = "bray")

UE.pair <- pairwise.adonis(UE.asin, UE$Fraction_size_um, sim.method = "bray")
DCM.pair <- pairwise.adonis(DCM.asin, DCM$Fraction_size_um, sim.method = "bray")
MESO.pair <- pairwise.adonis(MESO.asin, MESO$Fraction_size_um, sim.method = "bray")

#Run homogeneity of multivariate dispersions (how different are samples within the same group? Are they really similar aka tight clustering or can they be kinda different)
all.bc <- vegdist(ASV_asin,method="bray")
all.beta <- betadisper(all.bc, group = mydata$depth.bin)
all.beta2 <- betadisper(all.bc, group = mydata$Fraction_size_um)
