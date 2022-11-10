#! /usr/bin/env Rscript

# Gabriel Hoffman
# April 6, 2021
#
# Plot disk usage based on showquota
# Must be run on minerva login node

# /hpc/users/hoffmg01/build2/hpc_disk_usage/plot_showquota.R --project psychgen

library(getopt)
spec = matrix(c(
      'project', 'p', 1, "character"
    ), byrow=TRUE, ncol=4)
opt = getopt(spec)

# projects = "CommonMind epigenAD roussp01a va-biobank psychencode epigenBD roussp01b psychAD psychgen psychgen2"

# command in bash
cmd = paste0("for proj in ", opt$project, "; do showquota -p $proj arion; done | grep -v Usage | grep -v 'Last report update' | awk '{print $1, $2, $3, $4}' | sed 's/psychgenadmin/psychgenadmin /g' > file_size.tmp")

system( cmd )

# R code
########

suppressPackageStartupMessages({
library(ggplot2)
library(tidyverse)
library(data.table)
library(fs)
library(R.utils)
library(stringr)
library(lubridate)
})

# get date of report
dates = system("showquota -p CommonMind arion  | head -n 1 | cut -f4 -d' '", intern=TRUE)
date_main = ymd_hm(dates)


# read size file
df = fread('file_size.tmp')
colnames(df) = c("Project", "User", "Size", "Files")

# get human name for each user
df$Name = sapply(df$User, function(id){
	cmd = paste("finger", id, "| grep 'Name:' | awk '{$1=$2=$3=\"\"; print }'")

	trimws(system(cmd, intern=TRUE))
})

# convert sizes to numeric
df$Project = gsub("arion_projects_", '', df$Project)
df$Size = sapply(df$Size, fs_bytes) %>% 
	hsize(units="GB") %>% 
	gsub(pattern=" GB", replacement="") %>% 
	as.numeric
df = df[!is.na(Size),]

df_sort = df[,data.frame(TotalSize=sum(Size)), by="Name"]
df_sort = df_sort[order(TotalSize, decreasing=FALSE),]

df$Name = factor(df$Name, df_sort$Name)

# Big users
df_plot = df[Name %in% df_sort[TotalSize >= 1000, Name],]

main = paste0("Big users (", date_main, ')')
fig1 = ggplot(df_plot, aes(Name, Size/1024, fill=Project)) + geom_bar(stat="identity") + coord_flip() + scale_y_continuous(expand=c(.01, .01)) + theme_bw() + ylab("Size (Tb)") + xlab("Human name") + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=2)

# Small users
df_plot = df[Name %in% df_sort[(TotalSize < 1000) & TotalSize > 1, Name],]

main = paste0("Small users (", date_main, ')')

fig2 = ggplot(df_plot, aes(Name, Size, fill=Project)) + geom_bar(stat="identity") + coord_flip() + scale_y_continuous(expand=c(.01, .01)) + theme_bw() + ylab("Size (Gb)") + xlab("Human name") + ggtitle(main) + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=2)



# write to file
file = paste0("/hpc/users/hoffmg01/www/DiskUsage/DiskUsage_", opt$project, "_", year(date_main), "-", month(date_main), "-", day(date_main) ,".pdf")
pdf(file)
fig1
fig2
dev.off()

# https://hoffmg01.u.hpc.mssm.edu/DiskUsage/


