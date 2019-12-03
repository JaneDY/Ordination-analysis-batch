#!/usr/bin/env python
#coding: utf8
    
__author__ = 'yun.ding'
        
	    
def pcoaGroup(inType, input, infile, output, ingroup, method, pc, ps, ls, width, height):
    Rcmd = '''
library(vegan)
library(ggpubr)
\n'''
    if inType == 'AF':
        Rcmd += '''
data <- read.table("'''+input+'''",sep="\\t",head=T,comment.char = "",check.names = FALSE, row.names = 1)
dataT <- t(data)
dist <- vegdist(dataT, method="'''+method+'''")
dist <- as.matrix(dist)
\n'''   
    if inType == 'DF':
        Rcmd += '''
dist <- read.table("'''+input+'''",sep="\\t",head=T,comment.char = "",check.names = FALSE, row.names = 1)
\n'''
    Rcmd += ''' 
sd <- read.table("'''+ingroup+'''",head=T,sep="\\t",comment.char = "",check.names = FALSE)
rownames(sd) <- as.character(sd[,1])
sd[,2] <- as.character(sd[,2])

dist <- dist[as.character(sd[,1]),][,as.character(sd[,1])]
\n'''
    if inType == 'AF':
        Rcmd += '''
write.table(dist, "'''+output+'''/newdist.'''+infile+'''_with_'''+method+'''.xls", sep="\\t", col.names=NA)
\n'''
    if inType == 'DF':
        Rcmd += '''
write.table(dist, "'''+output+'''/newdist.'''+infile+'''", sep="\\t", col.names=NA)
\n'''
    Rcmd += '''
pc_num <- as.numeric(unlist(strsplit("'''+pc+'''","-")))
pc_x <- pc_num[1]
pc_y <- pc_num[2]

pca <- cmdscale(dist, k=3, eig=TRUE)
pc12 <- pca$points[,pc_num]
pc <- round(pca$eig/sum(pca$eig)*100,digits = 2)

pc12 <- as.data.frame(pc12)
colnames(pc12) <- c("pc_x","pc_y")
pc12['sample'] <- rownames(pc12)
colnames(sd) <- c("sample","group")
pc12 <- merge(pc12,sd,by="sample")

mex <- 0.2*abs(max(pc12$pc_x)-min(pc12$pc_x)) 
mey <- 0.2*abs(max(pc12$pc_y)-min(pc12$pc_y))

ave_x <- tapply(pc12$pc_x, pc12$group, mean)
ave_y <- tapply(pc12$pc_y, pc12$group, mean)
#sd_x <- tapply(pc12$pc_x, pc12$group, sd)
#sd_y <- tapply(pc12$pc_y, pc12$group, sd)
sd_x <- tapply(pc12$pc_x, pc12$group, sd)/sqrt(table(pc12$group))
sd_y <- tapply(pc12$pc_y, pc12$group, sd)/sqrt(table(pc12$group))


gpc12 <- as.data.frame(cbind(ave_x,ave_y,sd_x,sd_y))
gpc12$group <- row.names(gpc12)

mytheme <- theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
                 text = element_text(size=18),
                 plot.margin = unit(c(2, 2, 2, 2), "lines"))

p <- ggscatter(gpc12, x = "ave_x", y = "ave_y", 
               color = "group", size = '''+ps+''', palette = "npg")

p <- p + geom_errorbar(aes(ymin=ave_y-sd_y, ymax=ave_y+sd_y, color=group), width=mex*0.1, size=1.2) + 
  geom_errorbarh(aes(xmin=ave_x-sd_x, xmax=ave_x+sd_x, color=group), height=mey*0.1, size=1.2)

p <- p + ylab(paste0("PC",pc_y,"(",round(pc[pc_y],2),"%",")")) + 
  xlab(paste0("PC",pc_x,"(",round(pc[pc_x],2),"%",")"))
p <- p + geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash") + 
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")
p <- p + theme(legend.position = "none")

p <- p + geom_text(aes(label = gpc12$group), size='''+ls+''', colour = '#000000', vjust=0.5, hjust=0.5)

p <- p + mytheme

name <- paste0("'''+output+'''/pcoa_by_group_PC",pc_x,"_PC",pc_y,".pdf")
ggsave(name, p, width = '''+width+''', height = '''+height+''', limitsize = F )
\n'''
    Rcmd += '''
write.table(pca$points, "'''+output+'''/pcoa_by_group_sites.xls", sep="\\t", col.names=NA)
write.table(pca$eig/sum(pca$eig), "'''+output+'''/pcoa_by_group_importance.xls", sep="\\t")
\n'''
    open(output+'/cmd.pcoa.group.r', 'w').write(Rcmd)
    os.system('Rscript '+output+'/cmd.pcoa.group.r')
    return 'Done PCoA type: group'

def pcoaBox(inType, input, output, ingroup, method, pc, boxType, plabel, ell):
    Rcmd = '''
library(vegan)
#library(philentropy)
library(ggpubr)
library(reshape2)
\n'''
    if inType == 'AF':
        Rcmd += '''
data <- read.table("'''+input+'''",sep="\\t",head=T,comment.char = "",check.names = FALSE, row.names = 1)
dataT <- t(data)
dist <- vegdist(dataT, method="'''+method+'''")
dist <- as.matrix(dist)
\n'''
    if inType == 'DF':
        Rcmd += '''
dist <- read.table("'''+input+'''",sep="\\t",head=T,comment.char = "",check.names = FALSE, row.names = 1)
\n'''
    Rcmd += '''
sd <- read.table("'''+ingroup+'''",head=T,sep="\\t",comment.char = "",check.names = FALSE)
rownames(sd) <- as.character(sd[,1])
sd[,2] <- as.character(sd[,2])

dist <- dist[as.character(sd[,1]),][,as.character(sd[,1])]

pc_num <- as.numeric(unlist(strsplit("'''+pc+'''","-")))
pc_x <- pc_num[1]
pc_y <- pc_num[2]

pca <- cmdscale(dist, k=3, eig=TRUE)
pc12 <- pca$points[,pc_num]
pc <- round(pca$eig/sum(pca$eig)*100,digits = 2)

pc12 <- as.data.frame(pc12)
colnames(pc12) <- c("pc_x","pc_y")
pc12['sample'] <- rownames(pc12)
colnames(sd) <- c("sample","group")
pc12 <- merge(pc12,sd,by="sample")

mex <- 0.2*abs(max(pc12$pc_x)-min(pc12$pc_x)) 
mey <- 0.2*abs(max(pc12$pc_y)-min(pc12$pc_y))

mytheme <- theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 text = element_text(size=18))
                 #plot.margin = unit(c(1, 0.5, 1, 1), "lines"))

p <- ggscatter(pc12, x = "pc_x", y = "pc_y", 
               color = "group", shape = "group", palette = "npg", size=3)

p <- p + ylab(paste0("PC",pc_y,"(",round(pc[pc_y],2),"%",")"))
p <- p + xlab(paste0("PC",pc_x,"(",round(pc[pc_x],2),"%",")"))
p <- p + geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash")
p <- p + geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")
p <- p + theme(legend.position = "none")
\n'''
    if plabel == 'T':
        Rcmd += '''
p <- p + geom_text(aes(label = pc12$sample), size=3, colour = '#595959', vjust=1.5, hjust=-0.2)
\n'''
    if ell == 'T':
        Rcmd += '''
p <- p + stat_ellipse(aes(x = pc_x, y = pc_y, fill = group), geom = "polygon", alpha = 0.3, level = 0.95, segments = 100) + 
  stat_ellipse(aes(x = pc_x, y = pc_y, color = group), linetype = 2, level = 0.99)
\n'''
    Rcmd += '''
p <- p + mytheme

#dist1 <- dist
#dist1[lower.tri(as.data.frame(dist1))] <- NA
#dist1$sample <- rownames(dist1)
#dist1 <- melt(dist1)
#dist1$variable <- as.character(dist1$variable)
#dist1 <- dist1[dist1$sample!=dist1$variable,]
#dist1 <- na.omit(dist1)
#dist1 <- merge(dist1,sd,by='sample',all= F)
#dist1 <- merge(dist1,sd,by.x='variable',by.y='sample',all = F)
#dist1$group <- ifelse(dist1$group.x==dist1$group.y,paste(dist1$group.x),paste('Inter'))
#intra <- dist1[!dist1$group %in% "Inter",]
#cp <- combn(unique(intra$group),2)
		  
cp <- combn(unique(pc12$group),2)
comp <- list()
for(i in 1:ncol(cp)){
          comp[[i]] <- cp[,i]
}
\n'''
    Rcmd += '''
pr <- ggboxplot(pc12, x="group", y="pc_y", fill = "group", add = "jitter", palette = "npg") + 
   stat_compare_means(comparisons = comp, label = "p.signif")
pr <- pr + theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
                 axis.text.y = element_blank(),
                 axis.ticks.y= element_blank(),
                 axis.title.y = element_blank(),
                 legend.position = "none",
                 axis.text.x = element_text(size = 15),
                 axis.title.x = element_text(size = 20))
\n'''
    if boxType == 'right':
        Rcmd += '''
p <- p + theme(plot.margin = unit(c(2, 0.5, 1.9, 1.5), "lines"))
pr <- pr + theme(plot.margin = unit(c(2, 2.5, 0.75, 0), "lines"))
'''
        if inType == 'AF':
            Rcmd += '''
name <- paste0("'''+output+'''/pcoa_with_boxplot_PC",pc_y,".pdf")
\n'''
        else:
            Rcmd += '''
name <- paste0("'''+output+'''/pcoa_with_boxplot_PC",pc_y,".pdf")
\n'''
        Rcmd += '''
pdf(name, width = 11, height = 8)
ggarrange(p, pr, ncol = 2, nrow = 1, widths = c(2,1), heights = c(1,1))#, align = "h")
dev.off()
\n'''
    Rcmd += '''
pt <- ggboxplot(pc12, x="group", y="pc_x", fill = "group", add = "jitter", palette = "npg") + coord_flip() + 
   stat_compare_means(comparisons = comp, label = "p.signif")
pt <- pt + theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "none",
                 axis.text.y = element_text(size = 15, angle = 90),
                 axis.title.y = element_text(size = 20))
\n'''
    if boxType == 'top':
        Rcmd += '''
p <- p + theme(plot.margin = unit(c(0, 2.5, 1.5, 1.5), "lines"))
pt <- pt + theme(plot.margin = unit(c(2, 2.5, 0.5, 0.5), "lines"))
'''
        if inType == 'AF':
            Rcmd += '''
name <- paste0("'''+output+'''/pcoa_with_boxplot_PC",pc_y,".pdf")
\n'''
        else:
            Rcmd += '''
name <- paste0("'''+output+'''/pcoa_with_boxplot_PC",pc_x,".pdf")
\n'''
        Rcmd += '''
pdf(name, width = 8, height = 11)
ggarrange(pt, p, ncol = 1, nrow = 2, widths = c(1,1), heights = c(1,2))#, align = "v")
dev.off()
\n'''
    if boxType == 'both':
        Rcmd += '''
p0 <- ggplot() + theme(panel.background = element_blank(),
                       plot.margin = unit(c(0, 0, 0, 0), "lines"))
p <- p + theme(plot.margin = unit(c(0, 0, 2, 2), "lines"))#plot.margin = unit(c(0, 0, 1.5, 1.5), "lines"))
pr <- pr + theme(plot.margin = unit(c(0, 2, 0.75, 0.5), "lines"))#plot.margin = unit(c(0, 2.5, 1.5, 0), "lines"))
pt <- pt + theme(plot.margin = unit(c(2, 0, 0.5, 0.8), "lines"))#plot.margin = unit(c(2, 0, 0, 1.5), "lines"))
'''
        if inType == 'AF':
            Rcmd += '''
name <- paste0("'''+output+'''/pcoa_with_boxplot_PC",pc_x,"_PC",pc_y,".pdf")
\n'''
        else:
            Rcmd += '''
name <- paste0("'''+output+'''/pcoa_with_boxplot_PC",pc_x,"_PC",pc_y,".pdf")
\n'''
        Rcmd += '''
pdf(name, width = 10, height = 10)
ggarrange(pt, p0, p, pr, ncol = 2, nrow = 2, widths = c(2,1), heights = c(1,2))#, align = "hv")
dev.off()
\n'''
    Rcmd += '''
write.table(pca$points, "'''+output+'''/pcoa_with_boxplot_sites.xls", sep="\\t", col.names=NA)
write.table(pca$eig/sum(pca$eig), "'''+output+'''/pcoa_with_boxplot_importance.xls", sep="\\t")
\n'''
    open(output+'/cmd.pcoa.box.r', 'w').write(Rcmd)
    os.system('Rscript '+output+'/cmd.pcoa.box.r')
    return 'Done PCoA type: box'

def pcoaDiv(inType, input, infile, output, ingroup, method, adiv, index, pc, plabel, width, height):
    Rcmd = '''
library(vegan)
library(ggpubr)
library(ggsci)
\n'''
    if inType == 'AF':
        Rcmd += '''
da <- read.table("'''+input+'''",sep="\\t",head=T,check.names = FALSE)
rownames(da) <- as.character(da[,1])
da <- da[,-1]
da <- t(da)
dist <- vegdist(da, method="'''+method+'''")
dist <- as.matrix(dist)
\n'''
    if inType == 'DF':
        Rcmd += '''
dist <- read.table("'''+input+'''",sep="\\t",head=T,comment.char = "",check.names = FALSE, row.names = 1)
\n'''
    Rcmd += '''
indiv <- read.table("'''+adiv+'''",head=T,sep="\\t",comment.char = "",check.names = FALSE)
rownames(indiv) <- as.character(indiv[,1])
div <- indiv[,c(1,which(colnames(indiv)=="'''+index+'''"))]
\n'''
    if ingroup:
        Rcmd += '''
sd <- read.table("'''+ingroup+'''",head=T,sep="\\t",comment.char = "",check.names = FALSE)
rownames(sd) <- as.character(sd[,1])
sd[,2] <- as.character(sd[,2])
dist <- dist[as.character(sd[,1]),][,as.character(sd[,1])]
div <- div[as.character(sd[,1]),]
\n'''
    if inType == 'AF':
        Rcmd += '''
write.table(dist, "'''+output+'''/newdist.'''+infile+'''", sep="\\t", col.names=NA)
\n'''
    Rcmd += '''
pc_num <- as.numeric(unlist(strsplit("1-2","-")))
pc_x <- pc_num[1]
pc_y <- pc_num[2]

pca <- cmdscale(dist, k=3, eig=TRUE)
pc12 <- pca$points[,pc_num]
pc <- round(pca$eig/sum(pca$eig)*100,digits = 2)

pc12 <- as.data.frame(pc12)
colnames(pc12) <- c("pc_x","pc_y")
pc12['sample'] <- rownames(pc12)
\n'''	
    if ingroup:
        Rcmd += '''
colnames(sd) <- c("sample","group")
pc12 <- merge(pc12,sd,by="sample")

ave_x <- tapply(pc12$pc_x, pc12$group, mean)
ave_y <- tapply(pc12$pc_y, pc12$group, mean)

gpc12 <- as.data.frame(cbind(ave_x,ave_y))
gpc12$group <- row.names(gpc12)

merge_pc12 <- merge(pc12, gpc12, by="group")
pc12 <- merge_pc12
\n'''
    Rcmd += '''
colnames(div) <- c("sample","index")
pc12 <- merge(pc12,div,by="sample")

dmin <- min(pc12$index)
dmax <- max(pc12$index)

mex <- 0.2*abs(max(pc12$pc_x)-min(pc12$pc_x)) 
mey <- 0.2*abs(max(pc12$pc_y)-min(pc12$pc_y))

mytheme <- theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
                 axis.text.x = element_text(color = "black", size = 12),
                 axis.text.y = element_text(color = "black", size = 12),
                 text = element_text(size=18),
                 panel.background = element_blank(), 
                 plot.margin = unit(c(2.3, 4, 8, 3), "lines"))

mypch <-c(21:25,7:14)

p <- ggplot(pc12, aes(x=pc_x, y=pc_y, fill=index)) + geom_point(shape=21, size=4, stroke=1) + 
  scale_color_gsea() + 
  scale_fill_gsea(guide=guide_colorbar(barwidth = 13, direction = "horizontal", 
  title = "'''+index+''' diversity", title.position = "top"))

p <- p + guides(color=FALSE) + 
  theme(legend.title.align = 0.5, legend.position = c(0.65,0), legend.justification = c(1,1.6))
\n'''
    if ingroup:
	Rcmd += '''
p <- ggplot(pc12, aes(x=pc_x, y=pc_y, fill=index, color=index, shape=group)) + geom_point(size=4, stroke=0.1) + 
  scale_shape_manual(values = mypch[1:length(table(pc12$group))], guide=guide_legend(title = "group")) + scale_color_gsea() +
  scale_fill_gsea(guide=guide_colorbar(title = "'''+index+'''", title.position = "top"))

p <- p + guides(color=FALSE)
\n'''
    Rcmd += '''
p <- p + ylab(paste0("PC",pc_y,"(",round(pc[pc_y],2),"%",")")) + 
  xlab(paste0("PC",pc_x,"(",round(pc[pc_x],2),"%",")"))
p <- p + geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash") + 
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")
\n'''
    if plabel == 'T':
        Rcmd += '''
p <- p + geom_text(aes(label = pc12$sample), size=3, colour = '#595959', vjust=1.5, hjust=-0.2)
\n'''
#    if ingroup:
#        Rcmd += '''
#set <- unique(pc12$group)
#for(i in set){
#      p <- p + stat_ellipse(data = pc12[pc12$group==i,])
#}
#
#p <- p + geom_label(aes(x=ave_x, y=ave_y, label=group), fill="white", color="black", size=5)
#\n'''
    Rcmd += '''
p <- p + mytheme

name <- paste0("'''+output+'''/pcoa_with_'''+index+'''_PC", pc_x, "-PC", pc_y, ".pdf")
ggsave(name, p, width = '''+width+''', height = '''+height+''', limitsize = F)
\n'''
    Rcmd += '''
write.table(pca$points, "'''+output+'''/pcoa_with_adiv_sites.xls", sep="\\t", col.names=NA)
write.table(pc, "'''+output+'''/pcoa_with_adiv_importance.xls", sep="\\t")
\n'''
    open(output+'/cmd.pcoa.div.r', 'w').write(Rcmd)
    os.system('Rscript '+output+'/cmd.pcoa.div.r')
    return 'Done PCoA type: div'

def pcoaAbu(inType, input, infile, output, ingroup, method, abu, taxList, topTaxon, plabel, pc, width, height):
    Rcmd = '''
library(vegan)
library(ggpubr)
library(ggsci)
\n'''
    if inType == 'AF':
        Rcmd += '''
da <- read.table("'''+input+'''",sep="\\t",head=T,check.names = FALSE)
rownames(da) <- as.character(da[,1])
da <- da[,-1]
da <- t(da)
dist <- vegdist(da, method="'''+method+'''")
dist <- as.matrix(dist)
\n'''
    if inType == 'DF':
        Rcmd += '''
dist <- read.table("'''+input+'''",sep="\\t",head=T,comment.char = "",check.names = FALSE, row.names = 1)
\n'''
    Rcmd += '''
inabu <- read.table("'''+abu+'''",head=T,sep="\\t",comment.char = "",check.names = FALSE)
rownames(inabu) <- as.character(inabu[,1])
inabu <- inabu[,-1]
\n'''
    if ingroup:
        Rcmd += '''
sd <- read.table("'''+ingroup+'''",head=T,sep="\\t",comment.char = "",check.names = FALSE)
rownames(sd) <- as.character(sd[,1])
sd[,2] <- as.character(sd[,2])
dist <- dist[as.character(sd[,1]),][,as.character(sd[,1])]
inabu <- inabu[,as.character(sd[,1])]
\n'''
    if inType == 'AF':
        Rcmd += '''
write.table(dist, "'''+output+'''/newdist.'''+infile+'''", sep="\\t", col.names=NA)
\n'''
    Rcmd += '''
rowsum <- sapply(1:nrow(inabu),function(x) sum(inabu[x,]))
inabu <- inabu[order(rowsum,decreasing=TRUE),]
dat <- sapply(1:ncol(inabu),function(x) inabu[,x]/sum(inabu[,x])) #求相对丰度
colnames(dat) <- colnames(inabu)
rownames(dat) <- rownames(inabu)
\n'''
    if topTaxon:
        Rcmd += '''
topn <- 50
topTax <- dat[1:topn,]
\n'''
    elif taxList:
	    Rcmd += '''
tlist <- read.table("'''+taxList+'''", head=F,check.names = FALSE)
topTax <- dat[which(rownames(dat) %in% tlist[,1]),]
\n'''
    Rcmd += '''
abu <- as.data.frame(apply(topTax,2,sum))
colnames(abu) <- "abu"
abu$sample <- row.names(abu)

amin <- min(abu$abu)
amax <- max(abu$abu)
abuSeq <- seq(amin, amax, length.out = 3)
abuLabel <- sapply(abuSeq, function(x) paste0(round(x,2)*100,"%"))

pc_num <- as.numeric(unlist(strsplit("1-2","-")))
pc_x <- pc_num[1]
pc_y <- pc_num[2]

pca <- cmdscale(dist, k=3, eig=TRUE)
pc12 <- pca$points[,pc_num]
pc <- round(pca$eig/sum(pca$eig)*100,digits = 2)

pc12 <- as.data.frame(pc12)
colnames(pc12) <- c("pc_x","pc_y")
pc12['sample'] <- rownames(pc12)

pc12 <- merge(pc12,abu,by="sample")
\n'''
    if ingroup:
        Rcmd += '''
colnames(sd) <- c("sample","group")
pc12 <- merge(pc12,sd,by="sample")
\n'''
    Rcmd += '''
mex <- 0.2*abs(max(pc12$pc_x)-min(pc12$pc_x)) 
mey <- 0.2*abs(max(pc12$pc_y)-min(pc12$pc_y))

mytheme <- theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
                 axis.text.x = element_text(color = "black", size = 12),
                 axis.text.y = element_text(color = "black", size = 12),
                 text = element_text(size=18),
                 panel.background = element_blank(), 
                 plot.margin = unit(c(5, 4, 6, 2), "lines"))
\n'''
    if ingroup:
        Rcmd += '''
p <- ggscatter(pc12, x = "pc_x", y = "pc_y", palette = "npg") + 
  geom_point(aes(color=group,fill=group,size=abu), shape=21, stroke=2)
\n'''
    else:
        Rcmd += '''
p <- ggscatter(pc12, x = "pc_x", y = "pc_y") + geom_point(aes(size=abu), color="black", fill="black", shape=21, stroke=2)
\n'''
    Rcmd += '''
p <- p + ylab(paste0("PC",pc_y,"(",round(pc[pc_y],2),"%",")")) + 
  xlab(paste0("PC",pc_x,"(",round(pc[pc_x],2),"%",")"))
p <- p + geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash") + 
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")

p <- p + scale_size_continuous(range = c(3,8), limits = c(amin,amax), breaks = abuSeq, labels = abuLabel) + 
  theme(legend.position = "right") + 
  labs(size=NULL, color="Group", fill="Group") + 
  guides(fill = guide_legend(override.aes = list(size=3)))
\n'''
    if plabel == 'T':
        Rcmd += '''
p <- p + geom_text(aes(label = pc12$sample), size=3, colour = '#595959', vjust=1.5, hjust=-0.2)
\n'''
    Rcmd += '''
p <- p + mytheme

name <- paste0("'''+output+'''/pcoa_with_abund_PC", pc_x, "-PC", pc_y, ".pdf")
ggsave(name, p, width = '''+width+''', height = '''+height+''', limitsize = F)
\n'''
    Rcmd += '''
write.table(pca$points, "'''+output+'''/pcoa_with_abund_sites.xls", sep="\\t", col.names=NA)
write.table(pc, "'''+output+'''/pcoa_with_abund_importance.xls", sep="\\t")
\n'''
    open(output+'/cmd.pcoa.abu.r', 'w').write(Rcmd)
    os.system('Rscript '+output+'/cmd.pcoa.abu.r')
    return 'Done PCoA type: abu'

def pcoaRichness(inType, input, output, ingroup, method, otu, pc, plabel, width, height):
    Rcmd = '''
library(vegan)
library(ggpubr)
library(ggsci)
\n'''
    if inType == 'AF':
        Rcmd += '''
da <- read.table("'''+input+'''",sep="\\t",head=T,check.names = FALSE)
rownames(da) <- as.character(da[,1])
da <- da[,-1]
da <- t(da)
dist <- vegdist(da, method="'''+method+'''")
dist <- as.matrix(dist)
\n'''
    if inType == 'DF':
        Rcmd += '''
dist <- read.table("'''+input+'''",sep="\\t",head=T,comment.char = "",check.names = FALSE, row.names = 1)
\n'''
    Rcmd += '''
data <- read.table("'''+otu+'''",sep="\\t",head=T,comment.char = "",check.names = FALSE,row.name = 1)
\n'''
    if ingroup:
        Rcmd += '''
sd <- read.table("'''+ingroup+'''",head=T,sep="\\t",comment.char = "",check.names = FALSE)
rownames(sd) <- as.character(sd[,1])
sd[,2] <- as.character(sd[,2])
data <- data[,as.character(sd[,1])]
\n'''
    Rcmd += '''
# count otu number of each sample
number <- as.data.frame(apply(data, 2, function(x) length(which(x!=0))))
colnames(number) <- "count"
number$sample <- row.names(number)

cmin <- min(number$count)
cmax <- max(number$count)
countSeq <- seq(cmin, cmax, length.out = 4)
countLabel <- sapply(countSeq, function(x) round(x,0))

da <- t(data)
da <- da[,apply(da,2,function(x)any(x>0))]      

pc_num <- as.numeric(unlist(strsplit("1-2","-")))
pc_x <- pc_num[1]
pc_y <- pc_num[2]

pca <- cmdscale(dist, k=3, eig=TRUE)
pc12 <- pca$points[,pc_num]
pc <- round(pca$eig/sum(pca$eig)*100,digits = 2)

pc12 <- as.data.frame(pc12)
colnames(pc12) <- c("pc_x","pc_y")
pc12['sample'] <- rownames(pc12)
\n'''
    if ingroup:
        Rcmd += '''
colnames(sd) <- c("sample","group")
pc12 <- merge(pc12,sd,by="sample")
\n'''
    Rcmd += '''
pc12 <- merge(pc12,number,by="sample")

mex <- 0.2*abs(max(pc12$pc_x)-min(pc12$pc_x)) 
mey <- 0.2*abs(max(pc12$pc_y)-min(pc12$pc_y))

mytheme <- theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
                 text = element_text(size=18),
                 plot.margin = unit(c(6,2,5,2), "lines"))
\n'''
    if ingroup:
        Rcmd += '''
p <- ggscatter(pc12, x = "pc_x", y = "pc_y", 
               color = "group", size = "count", palette = "npg")
\n'''
    else:
        Rcmd += '''
p <- ggscatter(pc12, x = "pc_x", y = "pc_y", size = "count")
\n'''
    Rcmd += '''
p <- p + ylab(paste0("PC",pc_y,"(",round(pc[pc_y],2),"%",")")) + 
  xlab(paste0("PC",pc_x,"(",round(pc[pc_x],2),"%",")"))
p <- p + geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash") + 
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")
p <- p + scale_x_continuous(limits = c(min(pc12$pc_x)-mex,max(pc12$pc_x)+mex))
p <- p + scale_y_continuous(limits = c(min(pc12$pc_y)-mey,max(pc12$pc_y)+mey))
\n'''
    if ingroup:
        Rcmd += '''
p <- p + scale_size_continuous(range = c(3,10)) + 
  guides(size = guide_legend(title = "N.OTUs")) + theme(legend.position = "right") + 
  guides(color = guide_legend(title = "Group", override.aes = list(size=5)))
\n'''
    else:
        Rcmd += '''
p <- p + scale_size_continuous(range = c(3,10)) + 
  guides(size = guide_legend(title = "N.OTUs")) + theme(legend.position = "right")
\n'''
    if plabel == 'T':
        Rcmd += '''
p <- p + geom_text(aes(label = pc12$sample), size=3, colour = '#595959', vjust=1.5, hjust=-0.2)
\n'''
    Rcmd += '''
p <- p + mytheme

name <- paste0("'''+output+'''/pcoa_with_richness_PC", pc_x, "-PC", pc_y, ".pdf")
ggsave(name, p, width = '''+width+''', height = '''+height+''', limitsize = F)
\n'''
    Rcmd += '''
write.table(pca$x, "'''+output+'''/pcoa_with_richness_sites.xls", sep="\\t", col.names=NA)
write.table(summary(pca)$importance[2,], "'''+output+'''/pcoa_with_richness_importance.xls", sep="\\t")
\n'''
    open(output+'/cmd.pcoa.richness.r', 'w').write(Rcmd)
    os.system('Rscript '+output+'/cmd.pcoa.richness.r')
    return 'Done PCoA type: richness'

def paras():
    parser = argparse.ArgumentParser(description='Do Principal Coordinates Analysis.', epilog='What can I help U? Please contact yun.ding@majorbio.com.')
    parser.add_argument('-type', required = True, choices = ['group','box','div','abu','richness'], help = 'Choose which type of pcoa you want to plot')
    parser.add_argument('-mode', choices = ['AF', 'DF'], default = 'DF', help = 'Input file type, default: DF (AF: abundance file, DF: distance matrix)')
    parser.add_argument('-i', required = True, help = 'Input file')
    parser.add_argument('-o', default = './', help = 'Output directory, default: current directory')
    parser.add_argument('-g', help = 'Input group file, format: sampleID groupID')
    # for type abu
    parser.add_argument('-ab', help = 'Input taxon abundance file which for type abu')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-l', '--list', default = '', help = 'Taxon list, whose abundance will be used to set point size, with -ab')
    group.add_argument('-t', '--top', default = '50', help = 'Top abundance taxon to set point size, with -ab, default: 50')
    # for type div
    parser.add_argument('-div', help = 'Input diversity index file, format: sampleID    index1  index2  ...')
    parser.add_argument('-index', help = 'Choose which diversity index column to use, with -div')
    # for type richness
    parser.add_argument('-otu', help = 'Input otu table for type richness')
    # plot paras
    parser.add_argument('-s', choices = ['T', 'F'], default = 'F', help = 'Whether to show labels of points, default: F')
    parser.add_argument('-e', choices = ['T', 'F'], default = 'T', help = 'Whether to draw confidence ellipses, default: T')
    parser.add_argument('-pc', choices = ['1-2', '1-3', '2-3'], default = '1-2', help = 'Choose which PC axis to plot, default: 1-2')
    parser.add_argument('-m', choices = ["manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis"], default = 'bray', help = 'default: bray [When input an abundance file, then choose one dissimilarity index to generate the distance matrix]')
    parser.add_argument('-box', choices = ['right','top','both'], default = 'right', help = 'Choose which box to plot, default: right')
    parser.add_argument('-ps', default = '15', help = 'Point size, default: 15')
    parser.add_argument('-ls', default = '6', help = 'Label size, default: 6')
    parser.add_argument('-ww', default = '10', help = 'width, default: 10')
    parser.add_argument('-hh', default = '10', help = 'height, default: 10')
    args = parser.parse_args()
    return args


import os
import argparse

if __name__ == '__main__':

    args = paras()

    type = args.type
    inType = args.mode
    input = os.path.abspath(args.i)
    infile = os.path.split(input)[-1]
    output = os.path.abspath(args.o)
    if args.g:
        ingroup = os.path.abspath(args.g)
    else:
        ingroup = None
    method = args.m

    pc = args.pc
    ps = args.ps
    ls = args.ls
    boxType = args.box
    plabel = args.s
    ell = args.e
    width = args.ww
    height = args.hh

    if os.path.exists(output):
        print 'Warn: Output directory already exists'
    else:
        os.mkdir(output)
    
    if inType == 'AF':
	if not args.m:
	    print 'Error: Please use -m to set one dissimilarity index!'
	    exit()

    if type == 'group':
        if not ingroup:
            exit('Please input group file')
        else:
            pcoaGroup(inType, input, infile, output, ingroup, method, pc, ps, ls, width, height)
    if type == 'box':
        if not ingroup:
            exit('Please input group file')
        else:
            pcoaBox(inType, input, output, ingroup, method, pc, boxType, plabel, ell)
    if type == 'div':
        if not args.div or not args.index:
            exit('Please input diversity index file and choose one index by -div and -index')
        else:
            adiv = os.path.abspath(args.div)
            index = args.index
            pcoaDiv(inType, input, infile, output, ingroup, method, adiv, index, pc, plabel, width, height)
    if type == 'abu':
        if not args.ab:
            exit('Please input taxon abundance file by -ab')
        abu = os.path.abspath(args.ab)
        taxList = os.path.abspath(args.list)
        topTaxon = args.top
        pcoaAbu(inType, input, infile, output, ingroup, method, abu, taxList, topTaxon, plabel, pc, width, height)
    if type == 'richness':
	if not args.otu:
	    exit('Please input otu table by -otu')
	else:
	    otu = os.path.abspath(args.otu)
            pcoaRichness(inType, input, output, ingroup, method, otu, pc, plabel, width, height)

