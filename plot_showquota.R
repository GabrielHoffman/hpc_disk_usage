#! /usr/bin/env Rscript

# Gabriel Hoffman
# April 6, 2021
#
# Plot disk usage based on showquota
# Must be run on minerva login node

# /hpc/users/hoffmg01/build2/hpc_disk_usage/plot_showquota.R

# command in bash
cmd = "for proj in CommonMind epigenAD roussp01a va-biobank psychencode epigenBD roussp01b psychAD psychgen psychgen2; do showquota -p $proj arion; done | grep -v Usage | grep -v 'Last report update' | awk '{print $1, $2, $3, $4}' | sed 's/psychgenadmin/psychgenadmin /g' > file_size.tmp"

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
dates = system("grep 'Last report update' usage.wsgi | awk '{print $4}' | sort -u", intern=TRUE)
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
file = '/hpc/users/hoffmg01/www/DiskUsage/DiskUsage_showquota.pdf'
pdf(file)
fig1
fig2
dev.off()




