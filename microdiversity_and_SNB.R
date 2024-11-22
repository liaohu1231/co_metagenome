library(ggplot2)
library(ggsci)
library(ggpubr)
mic=read.table("all.scaffold_nucl_diver.tsv",sep = '\t',header = T,row.names = 1)
metadata=read.csv("../Macrodiversity/metagenome_env_urban.csv",row.names = 1)
rownames(metadata)=metadata$ID
length=read.table("GSV_CVD_LVD_95-80.length",sep = '\t',row.names = 1)
mic$length=length[as.character(rownames(mic)),1]
mic=na.omit(mic)
mic=mic[mic$length>5000,]
mic=mic[,-157]
##overall
ol=as.data.frame(matrix(NA,nrow = ncol(mic),ncol = 1))
colnames(ol)="Microdiversity"
rownames(ol)=colnames(mic)
for(i in 1:ncol(mic)){
  vl=mic[,i]
  vl2=mic[vl>0,i]
  vl2=na.omit(vl2)
  x=c()
  if ( length(vl2) >= 100 ){
    for (sa in 1:1000) {
      tmp=sample(vl2,size = 100,replace = F)
      tmp=mean(tmp)
      x=append(x,tmp)
    }
    am=mean(x)
  }else{
    am=mean(vl2)
  }
  ol$Microdiversity[i]=am
}
ol$landuses=metadata[rownames(ol),2]

ggplot(ol,aes(Microdiversity,landuses))+
  geom_violin(aes(fill=landuses),alpha=0.4)+
  geom_boxplot(aes(fill=landuses),width=0.2,
               outlier.size = 0)+
  theme_bw()+
  scale_fill_aaas()+
  stat_compare_means(method = "t.test",method.args = list(alternative="greater"),
                     comparisons = list(c("Forest","Residential area"),
                                        c("Residential area","Farmland"),
                                        c("Forest","Farmland")))+
  labs(x="Average microdiversity")
ggsave("a_mi_lu.pdf",device = "pdf",width = 6,height = 3.3)

ol$latitude=metadata[rownames(ol),8]
ol$longitude=metadata[rownames(ol),7]
ggplot(ol,aes(x = longitude,
                y = Microdiversity,group=landuses))+
  theme_classic()+
  geom_point(aes(color=landuses,fill=landuses,shape=landuses),
             size=2,alpha=0.5)+
  stat_smooth(aes(color=landuses,fill=landuses),
              alpha=0.1)+
  scale_color_aaas()+
  scale_fill_aaas()+
  scale_shape_manual(values = c(21,22,23,24))+
  scale_x_continuous(breaks = seq(90,145,by=5))+
  scale_y_continuous(breaks = seq(0,0.02,by=0.002),
                     limits = c(0.002,0.016))
ggsave("longitude_losses.pdf",device = 'pdf',width = 5,height = 3)
p2=ggplot(ol,aes(x = latitude,y = Microdiversity))+
  theme_classic()+
  geom_point(aes(color=landuses,fill=landuses,
                 shape=landuses),
             size=2,alpha=0.5)+
  stat_smooth(aes(group=landuses,
                  color=landuses,fill=landuses),
              alpha=0.1)+
  scale_color_aaas()+
  scale_fill_aaas()+
  scale_shape_manual(values = list(21,22,23,24))+
  scale_x_continuous(breaks = seq(20,50,by=5))+
  scale_y_continuous(breaks = seq(0,0.02,by=0.002),
                     limits = c(0.002,0.016))
p2
ggsave("latitude_losses.pdf",plot = p2,device = 'pdf',width = 5,height = 3)

##mantel
library(vegan)
library(dplyr)
metadata$microdiversity=ol[rownames(metadata),1]
mantel(xdis = vegdist(x = metadata$microdiversity,method = 'bray'),
       ydis = vegdist(x = metadata$pH,method = 'bray'),permutations = 10000,
       method = 'spearman',parallel = 5)
lm(formula = Moisture~microdiversity,metadata)%>%summary()
lm(formula = pH~microdiversity,metadata)%>%summary()
lm(formula = TP~microdiversity,metadata)%>%summary()
lm(formula = BIO1~microdiversity,metadata)%>%summary()

ggplot(metadata,aes(microdiversity,Moisture))+
  geom_point(aes(fill=Landuse),shape=21,size=2,alpha=0.6)+
  theme_bw()+
  geom_smooth(method = 'lm')+
  scale_fill_aaas()
ggsave("Moisture_MI.pdf",device = "pdf",width = 4.5,height = 2.2)
ggplot(metadata,aes(microdiversity,pH))+
  geom_point(aes(fill=Landuse),shape=21,size=2,alpha=0.6)+
  theme_bw()+
  geom_smooth(method = 'lm')+
  scale_fill_aaas()
ggsave("pH_MI.pdf",device = "pdf",width = 4.5,height = 2.2)
ggplot(metadata,aes(microdiversity,BIO1))+
  geom_point(aes(fill=Landuse),shape=21,size=2,alpha=0.6)+
  theme_bw()+
  geom_smooth(method = 'lm')+
  scale_fill_aaas()
ggsave("BIO1_MI.pdf",device = "pdf",width = 4.5,height = 2.2)
##social niche width
rpkm=read.table("../Macrodiversity/merged_GSV_CVD_LVD.rpkm",sep='\t',row.names = 1,header = T)
rpkm$length=length[rownames(rpkm),1]
rpkm=na.omit(rpkm)
rpkm=rpkm[rpkm$length>5000,]
rpkm=rpkm[,-157]
sum=apply(rpkm,1,sum)
rpkm=rpkm[sum>0,]
write.table(rpkm,"../Macrodiversity/merged_GSV_CVD_LVD.rpkm2",quote = F,row.names = T,sep = '\t')

relative_abundance=function(d){
  dta_sum=apply(d,2,function(x){x/sum(x)})
}
rpkm=relative_abundance(rpkm)
rpkm=as.data.frame(t(rpkm))


filter=function(data){
  for(i in 1:nrow(data)){
    for (j in 1:ncol(data)) {
      print(j)
      if(data[i,j] < 0.0001){
        data[i,j]=0
      }
    }
  }
}
library(reshape2)
SNB=function(rpkm2){
  v_snb=matrix(NA,nrow = ncol(rpkm2),ncol = 1)
  colnames(v_snb)="SNB"
  rownames(v_snb)=colnames(rpkm2)
  for(i in 1:ncol(rpkm2)){
    tmp=rpkm2[rpkm2[,i]>0,]
    print(i)
    sp = cor(x = t(tmp),method = "spearman")
    sp=melt(sp)
    sp$value2=0.5-sp$value/2
    r=sum(sp$value2)
    n=nrow(tmp)^2-nrow(tmp)
    SNB=r/n
    v_snb[i,]=SNB
  }
  v_snb=as.data.frame(v_snb)
  return(v_snb)
}
rpkm2=apply(rpkm, 2, function(x) ifelse(x < 0.0001, 0, x))
sum=apply(rpkm2,2,sum)
rpkm2=as.data.frame(rpkm2[,sum>0])
sum=apply(rpkm2,2,function(x){sum(x>0)})
rpkm2=as.data.frame(rpkm2[,sum>1])
v_snb=SNB(as.data.frame(rpkm2))
##average microdiversity
v_snb2=v_snb
v_snb2$mic=NA
for (i in 1:nrow(v_snb)){
  print(i)
  tmp=mic[rownames(mic)==rownames(v_snb2)[i],]
  if(nrow(tmp)>0){
    tmp[is.na(tmp)]=0
    tmp=as.data.frame(tmp[,tmp[1,]>0])
    s=list()
    for (n in 1:1000) {
      if(ncol(tmp)>5){
        t_m=mean(sample(t(tmp),5,replace = F))
        s=append(s,t_m)
      }else{
        t_m=mean(t(tmp)[,1])
        s=append(s,t_m)
      }
    }
    v_snb2$mic[i]=mean(as.data.frame(s)[,1])
  }
}
v_snb2=na.omit(v_snb2)
ggplot(v_snb2,aes(log10(mic),SNB))+
  geom_point(alpha=0.2,size=0.5,shape=21,fill="#6BF7E3")+
  theme_bw()
lm(SNB~log(mic),data = v_snb2)%>%summary()
cor.test(v_snb2$SNB,v_snb2$mic,method = "spearman")
ggsave("SNB_MI.pdf",device = "pdf",width = 4.5,height = 3.3)

tmp=strsplit(rownames(v_snb2),'_')
tmp=do.call(rbind,tmp)
v_snb2$sample=tmp[,1]
v_snb2$landuses=metadata[as.character(v_snb2$sample),2]
for(i in 1:nrow(v_snb2)){
  if( is.na(v_snb2$landuses[i])==TRUE ){
    tmp=strsplit(rownames(v_snb2)[i],'_')
    tmp=as.data.frame(do.call(rbind,tmp))
    if (tmp[1,2] == "resident"){
      v_snb2$landuses[i]="Residential area"
    }
    if(tmp[1,2] == "park"){
      v_snb2$landuses[i] = "Park"
    }
    if(tmp[1,2]=="forest"){
      v_snb2$landuses[i]="Forest"
    }
    if(tmp[1,2]=="farmland"){
      v_snb2$landuses[i]="Farmland"
    }
  }
}
v_snb2=na.omit(v_snb2)
ggplot(v_snb2,aes(log10(mic),SNB,color=landuses))+
  geom_point(alpha=0.4,size=1)+
  theme_bw()+
  scale_color_aaas()
ggsave("SNB_MI.pdf",device = "pdf",width = 5.5,height = 3.3)


ggplot(v_snb2,aes(landuses,SNB,color=landuses))+
  geom_boxplot(fill="grey90")+
  theme_bw()+
  scale_color_aaas()+
  stat_compare_means(comparisons = list(c("Forest","Farmland"),
                                        c("Forest","Park"),
                                        c("Farmland","Residential area")))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("SNB_box.pdf",device = "pdf",width = 5.5,height = 3.8)

v_snb[,2:157]=mic[rownames(v_snb),]
SNB_mi=matrix(NA,nrow = 156,ncol = 2)
colnames(SNB_mi)=c("R",'p')
rownames(SNB_mi)=colnames(mic)
SNB_mi=as.data.frame(SNB_mi)
for(i in 2:157){
  tmp=v_snb[,c(1,i)]
  tmp=tmp[tmp[,2]>0,]
  sp=cor.test(tmp[,1],tmp[,2],method = "pearson")
  SNB_mi$R[i-1]=sp$estimate
  SNB_mi$p[i-1]=sp$p.value
}
SNB_mi$landuses=metadata[rownames(SNB_mi),2]
SNB_mi$area=metadata[rownames(SNB_mi),3]
ggplot(SNB_mi,aes(landuses,R,color=landuses))+
  geom_boxplot()+
  theme_bw()+
  scale_color_aaas()+
  stat_compare_means(comparisons = list(c("Residential area","Park"),
                                        c("Residential area","Farmland")))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave("SNB_MI_spearman_box.pdf",device = "pdf",width = 5.5,height = 3.8)
ggplot(SNB_mi,aes(area,R,color=area))+
  geom_boxplot()+
  theme_bw()+
  scale_color_d3()+
  stat_compare_means(comparisons = list(c("Mid-temperate","Subtropical")))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

##sylph
sylph=read.table("sylph_global_soil_viruses_abu.tsv",sep='\t',header = 1)
sylph_ani=sylph[,c(1,15,13)]
library(reshape2)
sylph_ani=reshape2::acast(sylph_ani,formula = Contig_name~Sample_file)
sylph_ani_east=sylph_ani[,-grep("-",colnames(sylph_ani))]
sylph_ani_east=as.data.frame(sylph_ani_east)
sylph_ani_east$SNB=v_snb[rownames(sylph_ani_east),1]

SNB_ani=matrix(NA,nrow = 156,ncol = 2)
colnames(SNB_ani)=c("R",'p')
rownames(SNB_ani)=colnames(sylph_ani_east[,1:156])
SNB_ani=as.data.frame(SNB_ani)
for(i in 1:156){
  tmp=sylph_ani_east[,c(i,157)]
  tmp=na.omit(tmp)
  sp=cor.test(tmp[,1],tmp[,2],method = "pearson")
  SNB_ani$R[i]=sp$estimate
  SNB_ani$p[i]=sp$p.value
}

sylph_abu=sylph[,c(1,15,3)]
sylph_abu=reshape2::acast(sylph_abu,formula = Contig_name~Sample_file)
sylph_abu=as.data.frame(sylph_abu)
sylph_abu[is.na(sylph_abu)]=0
library(vegan)
distance <- as.matrix(vegdist(t(sylph_abu[,-grep("-",colnames(sylph_abu))]), method= "bray",na.rm = T))
pcoa_v <- cmdscale(distance, k = (nrow(t(sylph_abu[,-grep("-",colnames(sylph_abu))])) - 1), eig = TRUE)
pcoa_eig_v <- (pcoa_v$eig)[1:2] / sum(pcoa_v$eig)
sample_site_v <- data.frame ({pcoa_v$point})[1:3]
sample_site_v$names <- rownames(sample_site_v)
names(sample_site_v)[1:3] <- c('PCoA1', 'PCoA2','PCoA3')

metadata=read.csv("metadata2.txt",header = 1,sep='\t',row.names = 1)
sample_site_v$area=metadata[rownames(sample_site_v),5]
sample_site_v$types=metadata[rownames(sample_site_v),2]
library(ggplot2)
library(ggsci)
p <- ggplot(sample_site_v, aes(PCoA1, PCoA2)) +
  theme_bw() +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(aes(color=area,fill=area,shape=area),size = 2, alpha = 0.5) + 
  scale_fill_d3(palette = "category20") + 
  scale_color_d3(palette = "category20") + 
  scale_shape_manual(values = c(19,20,21,22,23,24,25))+
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig_v[1], 2), '%'), 
       y = paste('PCoA axis2: ', round(100 * pcoa_eig_v[2], 2), '%'))+
  #stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = area),
               #level = 0.95, show.legend = TRUE,inherit.aes = F)+
  #stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = site),level = 0.975, show.legend = TRUE,inherit.aes = F)+
  annotate('text', colour="#8766d")
p

p_value_1=adonis2(formula =distance~area,sample_site_v,permutations = 9999)
p_value_1

p_value_2=adonis2(formula = distance~types,sample_site_v,permutations = 9999)
p_value_2

p_value=anosim(t(sylph_abu),sample_site_v$area,permutations = 9999)
summary(p_value)