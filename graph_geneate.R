graph.generate <- function (file_name1, file_name2){
  text_directory <- c("file_name1")
  raw_text <- tolower(readChar(text_directory, nchars=10e6, useBytes = FALSE))
  raw_text <- gsub("[\t\r;]","",raw_text)
  raw_text <- gsub("\n"," ",raw_text)
  text <- gsub("\\[|\\]", "", raw_text)
  text <- removePunctuation(text, preserve_intra_word_dashes = TRUE)
  
  units <- length(strsplit(text, " ")[[1]])
  
  key_directory <- c("file_name2")
  raw_key <- tolower(readChar(key_directory, nchars=10e6, useBytes = FALSE))
  key <- gsub("[\t\r;]","", raw_key)
  key <- gsub("\n", " ",key)
  key <- gsub("\\[|\\]", "", key)
  key_list <- gsub("[[:punct:]]", "", key)
  key_list <- strsplit(key_list, " ")[[1]]
  
  out <- fcm(tokens(text), context = "window", window = 2, tri=F)
  n <- nrow(out)
  test <- as.matrix(cbind(out%*%rep(1,n),seq(n)))
  A <- as.matrix(out)
  D <- diag(rowSums(A))
  degree <- sum((as.matrix(out)!=0)%*%rep(1,n))
  truth <- unique(test[key_list[key_list %in% rownames(test)],2])
  y <- rownames(test) %in% key_list

  # get observed keywords from a third file or other resources
  #observed_keys <- 
  observed_keys <- paste(observed_keys, collapse = " ")
  
  # output: number of candidates, weighted adjacency matrix, degree matrix
  return(list("n"=n,"A"=A,"D"=D,"test"=test,"truth"=truth, "y" = y, "units"=units,"degree"=degree,
              key_list = key_list, observed_keys = observed_keys))
}