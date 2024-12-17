library(quanteda)
library(udpipe)
library(dplyr)
library(SnowballC)
library(tm)
library(quanteda.textstats)
library(caret)
library(MASS)
library(matrixcalc)

udmodel <- udpipe_download_model(language = "english-ewt")
udmodel <- udpipe_load_model(file = udmodel$file_model)

text_process <- function (dir1, dir2){
  directory <- dir1
  text <- tolower(readChar(directory, nchars=10e6, useBytes = FALSE))
  raw_text <- tolower(readChar(directory, nchars=10e6, useBytes = FALSE))
  raw_text<- gsub("[\t\r;]"," ",raw_text)
  raw_text<- gsub("\n"," ",raw_text)
  raw_text <- gsub("[/=]", " ", raw_text)
  raw_text <- removePunctuation(raw_text, preserve_intra_word_dashes = TRUE)
  x <- udpipe_annotate(udmodel,  raw_text)
  x <- as.data.frame(x)
  x <- x %>% dplyr::select(token, upos) %>% filter(upos %in% c("NOUN", "VERB", "ADJ"))
  tokenized_text <- x %>% mutate(stem = wordStem(x$token))
  tokenized_text <- tokenized_text$stem
  tokenized_text <- removePunctuation(tokenized_text)
  tokenized_text <- gsub("[[:punct:]]", "", tokenized_text)
  tokenized_text <- tokenized_text[tokenized_text != ""]
  temp <- textstat_frequency(dfm(tokens(tokenized_text)))
  tokenized_text <- tokenized_text[tokenized_text %in% temp$feature[which(temp$frequency > 2)]]
  tokenized_text <- paste(tokenized_text, collapse = " ")
  
  
  key_directory <- dir2
  raw_key <- tolower(readChar(key_directory, nchars=10e6, useBytes = FALSE))
  key <- gsub("[\t\r;]","", raw_key)
  key <- gsub("\n", " ",key)
  key <- gsub("\\[|\\]", "", key)
  x <- udpipe_annotate(udmodel,  key)
  x <- as.data.frame(x)
  x <- x %>% dplyr::select(token, upos) %>% filter(upos %in% c("NOUN", "VERB", "ADJ"))
  tokenized_key <- x %>% mutate(stem = wordStem(x$token))
  
  out <- fcm(tokens(tokenized_text, context = "window", window = 2,tri=F))
  n <- nrow(out)
  test <- as.matrix(cbind(out%*%rep(1,n),seq(n)))
  
  key_list <- unique(tokenized_key$stem)
  key_list <- key_list[key_list %in% rownames(test)]
  key_list <- paste(key_list, collapse = " ")
  key_list <- removePunctuation(key_list)
  key_list <- gsub("[[:punct:]]", "", key_list)
  
  len_key <- length(strsplit(key_list, " ")[[1]])
  
  return(list(doc = tokenized_text, keywords = key_list, len_key = len_key))
}
