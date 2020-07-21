# April 23, 2019
#
# Gabriel Hoffman
#
# Get disk usage for each user in groups: CommonMind, psychencode, psychgen
# Plot summary of results

# Get disk usage per user
# aggregate directly in awk

PROJECTS=('/sc/arion/projects/CommonMind' '/sc/arion/projects/psychencode' '/sc/hydra/projects/psychgen' '/sc/hydra/projects/roussp01a' '/sc/hydra/projects/roussp01b' '/sc/hydra/projects/epigenAD' '/sc/arion/projects/psychAD' '/sc/hydra/projects/va-biobank')

OUT=~/work/disk_usage_all.tsv
OUTBAM=~/work/disk_usage_BAM.tsv
OUTFQ=~/work/disk_usage_FQ.tsv
OUTRAWFQ=~/work/disk_usage_rawFQ.tsv
OUTRAWVCF=~/work/disk_usage_rawVCF.tsv

\rm -f $OUT $OUTBAM $OUTFQ $OUTRAWFQ $OUTRAWVCF

echo -e "Project\tUser\tSize" > $OUT
for GRP in ${PROJECTS[@]};
do
	# total size
	find $GRP -type f -printf "%u %s\n" | awk -v grp=$GRP -v OFS='\t' '{a[$1] += $2} END{for (i in a) print grp, i, a[i]}' >> $OUT

	# BAM
	find $GRP -name "*.bam" -type f -printf "%u %s\n" | awk -v grp=$GRP -v OFS='\t' '{a[$1] += $2} END{for (i in a) print grp, i, a[i]}' >> $OUTBAM

	# FASTQ
	find $GRP -name "*fastq.gz" -type f -printf "%u %s\n" | awk -v grp=$GRP -v OFS='\t' '{a[$1] += $2} END{for (i in a) print grp, i, a[i]}' >> $OUTFQ

	# raw FASTQ
	find $GRP -name "*fastq" -type f -printf "%u %s\n" | awk -v grp=$GRP -v OFS='\t' '{a[$1] += $2} END{for (i in a) print grp, i, a[i]}' >> $OUTRAWFQ

	# raw VCF
	find $GRP -name "*vcf" -type f -printf "%u %s\n" | awk -v grp=$GRP -v OFS='\t' '{a[$1] += $2} END{for (i in a) print grp, i, a[i]}' >> $OUTRAWVCF
done


# get user information
cat $OUT | grep -v "^Project" | cut -f2 | sort -u | parallel -P1 "finger {}" | grep Name | cut -f1,4 | sed 's/Login: //g' | sed 's/Name: //g' | tail -n +2  |  tr -s ' ' | tr '\t' ',' | sed 's/ ,/,/g' > ~/work/user_names.tsv





ll /sc/arion/projects/CommonMind/roussp01a/fetal_hic/hicpro/step1
ll /sc/arion/projects/CommonMind/roussp01a/iPSC_goate/step1/files/ips4
ll /sc/arion/projects/CommonMind/yixuan
ll /sc/arion/projects/CommonMind/INGELHEIM/atacseq/step1/files/942_N/star/94
ll /sc/arion/projects/CommonMind/INGELHEIM/rnaseq/step1/work/4a
ll /sc/hydra/projects/psychgen/sweden/swedex2/data/genome/out7-16
ll /sc/hydra/projects/roussp01a/bzeng/Capstone4_ROSMAP/trans-eQTL_detection/

/sc/hydra/projects/roussp01b/Georgios/CMC/fastq.tables/DLPFC.txt




# DOUG


You have 250,943 UNZIPPED VCFs

df = fread("dougs_file_usage.txt")
colnames(df) = c("User", "File", "Size")

df = df[order(df$Size, decreasing=TRUE),]
df$Size = df$Size/1e9
head(df)


find /sc/hydra/projects/psychgen -name "*vcf" -type f -printf "%u %p %s\n" | grep ruderd02 > dougs_file_usage.txt

/sc/hydra/projects/psychgen/sweden/swedexAutism/VCFs/header.Broad_families.combined.anno.vcf




##########################
# R code to plot results #
##########################

library(ggplot2)
library(gridExtra)
library(foreach)
library(data.table)

userNames = read.csv('~/work/user_names.tsv', stringsAsFactors=FALSE)
colnames(userNames) = c("login", "name")
userNames$login = trimws(userNames$login )

files = c(All = '~/work/disk_usage_all.tsv', BAM = '~/work/disk_usage_BAM.tsv', FASTQ = '~/work/disk_usage_FQ.tsv', RAW_FASTQ = '~/work/disk_usage_rawFQ.tsv', RAW_VCF='~/work/disk_usage_rawVCF.tsv')

df = lapply( 1:length(files), function(i){

	cmd = paste("cat ", files[i]," | awk '{if(NF==3) print $0}' | grep -v '^Project' > ~/work/disk_usage2.tsv")
	system( cmd )

	df = read.table('~/work/disk_usage2.tsv', stringsAsFactors=FALSE, header=FALSE)
	colnames(df) = c("Project", "User", "Size")
	df = data.table(df)
	df[,Total:=sum(Size), by="User"]
	df[,Metric := names(files)[i]]

	idx = match( df$User, userNames$login )
	df[,Name := userNames$name[idx]]
	df[,NameCombine := paste0(df$Name, ' (', df$User, ')')]

	df
})
df = do.call(rbind, df)



# keep individuals with > 10 Gb 
df_keep = df[,data.frame(keep = sum(Total)/1e9 > 100), by='Name']
df = df[Name %in% df_keep[keep==TRUE,Name],]


pdf("~/work/DiskUsage.pdf", height=16, width=16)
lapply( names(files), function(ID){
	ggplot(df[Metric==ID,], aes(reorder(NameCombine, Total), Size/1e12, fill=basename(Project))) + geom_bar(stat='identity') + theme_bw(16) + coord_flip() + ylab("Storage (Tb)") + ggtitle(paste(ID, "Storage:", date())) + theme(plot.title = element_text(hjust = 0.5)) + xlab("User")
})
lapply( names(files), function(ID){
	ggplot(df[Metric==ID,], aes(reorder(NameCombine, Total), Size/1e12, fill=basename(Project))) + geom_bar(stat='identity') + theme_bw(16) + coord_flip() + ylab("Storage (Tb)") + ggtitle(paste(ID, "Storage:", date())) + theme(plot.title = element_text(hjust = 0.5), legend.position="none") + xlab("User") + facet_wrap(~basename(Project), nrow=1, scales="free_x") 
})
dev.off()



# access plot at 
https://hoffmg01.u.hpc.mssm.edu//DiskUsage.pdf





# # df = foreach( file = files, .combine=rbind) %do% {

# 	title = gsub( "\\.txt", "", file)

# 	df = read.table(file, stringsAsFactors=FALSE, header=TRUE)
# 	df = data.frame(User = rownames(df), Storage = df[,'Size']/1e12, Project = title, stringsAsFactors=FALSE)
# # }	

# df = data.table(df)
# df[,Total:=sum(Storage), by="User"]

# idx = match( df$User, userNames$login )
# df$Name = userNames$name[idx]
# df$NameCombine = paste0(df$Name, ' (', df$User, ')')

# fig = ggplot(df, aes(reorder(NameCombine, Total), Storage, fill=Project)) + geom_bar(stat='identity') + theme_bw(16) + coord_flip() + ylab("Storage (Tb)") + ggtitle("Combined Storage") + theme(plot.title = element_text(hjust = 0.5))

# ggsave(file='~/work/disk_usage/DiskUsage.pdf', fig, height=12, width=10)


# # get active users
# cat user_logins.txt | cut -f1 -d' ' | sort -u | parallel -P1 "finger {}" | grep Name | cut -f1,4 | sed 's/Login: //g' | sed 's/Name: //g' > user_logins_full.txt

# df = read.table('~/work/disk_usage/user_logins_full.txt', sep='\t', stringsAsFactors=FALSE)
# colnames(df) = c("ID", "Name")
# df$ID = trimws(df$ID )
# df$Name = trimws(df$Name )
# df = df[,2:1]

# pdf('~/work/disk_usage/Logins.pdf', height=12, width=10)
# grid.table(df, row=NULL)
# dev.off()
