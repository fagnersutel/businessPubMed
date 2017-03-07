extract_pubMed_data <-
function (pubMed_query, batch_size = 1000, affi_regex_exclude = NULL) 
{
  ptm <- proc.time()
  my.idlist <- get_pubmed_ids(pubMed_query)
  record.num <- my.idlist$Count
  my.seq <- seq(1, as.numeric(my.idlist$Count), by = batch_size)
  pubmed.data <- lapply(my.seq, (function(ret.start) {
    batch.xml <- NULL
    message(paste("round #", which(my.seq == ret.start), 
                  " of ", length(my.seq), " ", sep = ""), appendLF = FALSE)
    while (is.null(batch.xml)) {
      batch.xml <- tryCatch({
        tmp.idlist <- get_pubmed_ids(pubMed_query)
        fetch_pubmed_data(tmp.idlist, retstart = ret.start, 
                          retmax = batch_size)
      }, error = function(e) {
        NULL
      })
    }
    record.list <- easyPubMed::articles_to_list(batch.xml)
    xtracted.data <- lapply(1:length(record.list), (function(i) {
      if (length(record.list) > 60) {
        custom.seq <- as.integer(seq(1, length(record.list), 
                                     length.out = 50))
        if (i %in% custom.seq) {
          message(".", appendLF = FALSE)
        }
      } else {
        message(".", appendLF = FALSE)
      }
      tmp.record <- tryCatch(article_to_df_dev(pubmedArticle = record.list[[i]], 
                                               autofill = TRUE, max_chars = 10), 
                             error = function(e) {
                               NULL
                             })
      if (!is.null(tmp.record)) {
        required.cols <- c("title", "year", "journal", 
                           "lastname", "firstname", "address", "email")
        out.record <- data.frame(matrix(NA, nrow = nrow(tmp.record), 
                                        ncol = length(required.cols)))
        colnames(out.record) <- required.cols
        match.cols <- colnames(tmp.record)[colnames(tmp.record) %in% 
                                             required.cols]
        out.record[, match.cols] <- tmp.record[, match.cols]
      } else {
        out.record <- NULL
      }
      out.record
    }))
    xtracted.data <- do.call(rbind, xtracted.data)
    message(" Filtering... ", appendLF = FALSE)
    xtracted.data <- xtracted.data[!is.na(xtracted.data$address), 
                                   ]
    if (!is.null(affi_regex_exclude)) {
      xtracted.data <- xtracted.data[regexpr(affi_regex_exclude, 
                                             xtracted.data$address,
                                             ignore.case = TRUE) < 0, ]
    }
    message("done!", appendLF = TRUE)
    xtracted.data
  }))
  stop.watch <- proc.time() - ptm
  pubmed.data <- do.call(rbind, pubmed.data)
  out.data <- list()
  out.data$params <- list()
  out.data$params$query_string <- pubMed_query
  out.data$params$pubMed_id_list <- my.idlist
  out.data$params$batch_size <- batch_size
  out.data$params$affi_regex_exclude <- affi_regex_exclude
  out.data$params$timing <- stop.watch
  out.data$data <- pubmed.data
  return(out.data)
}
