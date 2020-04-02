# Helper parameters and functions
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# plus minus symbol
plusminus <- "\u00b1"

#standard error
se <- function(x) {sqrt(var(x)/length(x))}

#default ggplot colors for two groups
d_colors <- c("#00BFC4", #lightblue
              "#F8766D") #red

#distinct colors
dist_colors <- c("#e6194b"
                 ,"#3cb44b"
                 ,"#ffe119"
                 ,"#0082c8"
                 ,"#f58231"
                 ,"#911eb4"
                 ,"#46f0f0"
                 ,"#f032e6"
                 ,"#d2f53c"
                 ,"#fabebe"
                 ,"#008080"
                 ,"#e6beff"
                 ,"#aa6e28"
                 ,"#fffac8"
                 ,"#800000"
                 ,"#aaffc3"
                 ,"#808000"
                 ,"#ffd8b1"
                 ,"#000080"
                 ,"#808080"
                 ,"#FFFFFF"
                 ,"#000000")

# default plot style
theme1 <- theme(text = element_text(family="Arial"),
                axis.text = element_text(size = 8),
                axis.title = element_text(size = 10),
                axis.text.x = element_text(margin = margin(10,0,10,0)),
                axis.text.y = element_text(margin = margin(0,10,0,10)),
                axis.title.x = element_text(margin = margin(0,0,10,0)),
                axis.title.y = element_text(margin = margin(0,0,0,10)))

# PLOS dimension rules
plos.link = "http://journals.plos.org/plosone/s/figures"
plos.imgtype = "tiff"
plos.resolution = 300
plos.filesize=10
plos.font = "Arial"
plos.fontsize = 10
plos.min.width = 6.68
plos.max.width = 19.05
plos.max.height = 22.23

plos_dims = function(width,height,scale = 1) {
  height_ratio = plos.max.height/height 
  width_ratio = plos.max.width/width
  use_ratio = 0
  if (height_ratio < width_ratio) {
    use_ratio = height_ratio
  } else {
    use_ratio = width_ratio
  }
  list(plos.height = height * use_ratio * scale,
       plos.width = width * use_ratio * scale)
}

# concatenates dplyr-grouped vectors into a 'list of vectors'
concat_grp <- function(column,isNum=F) {

  column = paste0(column,collapse=",")
  column = sapply(column,str_split,pattern=',')
  if (isNum) {column = lapply(column,as.numeric)}
  column
}

# collapses a 'list of vectors' into a vector of comma-separated values
unlist_row <- function(column,sep=",") {
  unlist(column) %>% paste(collapse=sep)
}

my_ttest <- function(pt,ctrl,show="pval") {
  
  # if patient and control values are all equal
  if (all(pt == ctrl)) {
    if (show == "higher") {
      "neither"
    } else {
      1
    }
  } else { # if patient and control values are not all equal
    res = t.test(x=pt,y=ctrl,paired=T,alternative = "two.sided")
    higher = ifelse(res$estimate > 0,"patient","control")
    
    if (show == "higher") {
      names(higher)=NULL
      higher
    } else {
      res$p.value
    }
    
  }
}

my_wsrtest <- function(pt,ctrl,show="pval") {
  
  # if patient and control values are all equal
  if (all(pt == ctrl)) {
    if (show == "higher") {
      "neither"
    } else {
      1
    }
  } else { # if patient and control values are not all equal
    res = coin::wilcoxsign_test(pt ~ ctrl,paired=T,zero.method="Pratt")
    higher = ifelse(res@statistic@teststatistic > 0,"patient","control")
    
    if (show == "higher") {
      names(higher)=NULL
      higher
    } else {
      res %>% pvalue()
    }
  }
}


which_higher <- function(comma_separated_num,
                         comma_separated_grp = c("patient","control")) {
  csm <- comma_separated_num
  csg <- comma_separated_grp
  csm <- str_split(csm,",") %>% unlist() %>% as.numeric()
  csg <- str_split(csg,",") %>% unlist()
  higher <- ifelse(all(csm[1]==csm),
                   "neutral",
                   csg[which.max(csm)])
  higher
}
# example
which_higher(comma_separated_num = c(1223,20),
             comma_separated_grp = c("patient","control"))


# if a regex is detected within a string, replaces the entire string with the replacement
replace_if <- function(str_vector, true_false,replacement) {
  str_vector[true_false] <- replacement
  str_vector
}

compare_counts <- function(count_df,count_column,grouping_column=NULL) {
  library(coin)
  # identify the counts column for comparison
  count_df <- count_df %>% mutate_("count" = count_column)
  
  # identify column to group the results (if not supplied)
  df_cols <- colnames(count_df)
  grouping_cols <- c(grouping_column,"custom","otu","species","genus","family","order","class","phyla")
  
  if (any(df_cols %in% grouping_cols)) {
    taxa_lvl <- df_cols[df_cols %in% grouping_cols][1]
  } else {
    stop("please specify the grouping column: i.e.'phyla','species','otu','gene','custom'")
  }


  # reshape dataframe for tests
  df_fortest <- count_df %>%
    arrange(stoneFormer,pair) %>%
    group_by_("stoneFormer",taxa_lvl) %>%
    summarize(pair = concat_grp(pair,isNum = T),
              count = concat_grp(count,isNum=T)) %>%
    spread(stoneFormer,count)

  # run ttest
  df_ttest <- df_fortest %>%
    mutate(test = "ttest") %>%
    mutate(pval = mapply(my_ttest,
                         pt=patient,
                         ctrl=control),
           higher = mapply(my_wsrtest,
                           pt=patient,
                           ctrl=control,show = "higher"))
  # run wsr test
  df_wsrtest <- df_fortest %>%
    mutate(test = "wilcoxon_signed_ranked") %>%
    mutate(pval = mapply(my_wsrtest,
                         pt=patient,
                         ctrl=control),
           higher = mapply(my_wsrtest,
                           pt=patient,
                           ctrl=control,show = "higher"))

  # combine test results into single df
  df_tests <- df_ttest %>%
    rbind(df_wsrtest) %>%
    select_(taxa_lvl,"test","pval","higher","everything()") %>%
    mutate(pair = sapply(pair,unlist_row),
           patient = sapply(patient,unlist_row),
           control = sapply(control,unlist_row)) %>%
    arrange(test,pval)

  #output
  df_tests
}

# helper function to remove NAs if there are other characters when pasting but defaulting back to NA if there are no other characters
# i.e. c(NA,"Hi","bob") -> "Hibob" #remove NA if there are other characters in the vector
# i.e. c(NA,NA) -> NA   #but default back to NA if no other characters
# i.e. c("Hi") -> "Hi"

paste_not_na <- function(vec,collapse=NULL) {
  if (all(is.na(vec))) {
    vec <- ""
    vec[vec ==""] <- NA
  } else {
    vec <- vec[!is.na(vec)]
    vec <- paste(vec,collapse = collapse)
  }
  vec
}


detect_str_in_vector <- function(str,str_vec,value=F,ignore_case=T) {
  detected <- str_detect(str,regex(str_vec,ignore_case = ignore_case))
  index <- which(detected)
  values <- str_vec[index]
  
  if (value) {
    values
  } else {
    length(index) > 0
  }
}
detect_str_in_vector("g__Oxalobacter",c("aple","oxalobacter","animalis"),value=F)
detect_str_in_vector("g__Oxalobacter",c("aple","oxalobacter","bobad"),value=T)
detect_str_in_vector("g__Oxalobacter",c("aple","oxalobacter","oxalo"),value=T)


count_gz_seqs <- function(gz_file,expected_filetype = ".fastq") {
  
  base_filename <- str_replace(gz_file,".gz$","")
  fastx_filename <- paste0(base_filename,
                           ifelse(str_detect(base_filename,expected_filetype),
                                  "",
                                  expected_filetype)
  )
  gunzip_cmd <- paste0("gunzip -c ",
                       gz_file," > ",
                       fastx_filename)
  system(gunzip_cmd)
  
  count_cmd <- paste("count_seqs.py -i", fastx_filename)
  output <- system(count_cmd,intern = T)
  output <- output %>% last() %>% str_extract("^[0-9]+") %>% as.numeric()
  
  system(paste("rm",fastx_filename))
  output
  
}
