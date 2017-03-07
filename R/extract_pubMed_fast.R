extract_pubMed_fast <-
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
    #
    all_titles   <- XML::xpathApply(batch.xml, "//ArticleTitle", saveXML)
    all_journals <- XML::xpathApply(batch.xml, "//Journal", saveXML)
    all_authors  <- XML::xpathApply(batch.xml, "//AuthorList", saveXML)
    #
    # length(all_titles)
    # length(all_journals)
    # length(all_authors)
    #
    if (length(all_titles) == length(all_journals) &
        length(all_journals) == length(all_authors))
    {
      data.out <- lapply(1:length(all_titles), (function(i){
        #
        if (length(all_titles)>100){
          if (i %in% as.integer(seq(1, length(all_titles), length.out = 50 )))
            message(".", appendLF = FALSE)
        }
        # Record specific info
        tmp.title <- custom_grep(all_titles[[i]], "ArticleTitle", format = "char")
        if (is.null(tmp.title)) { tmp.title <- NA }
        tmp.year  <- custom_grep(all_journals[[i]], "Year", format = "char")
        if (is.null(tmp.year)) { tmp.year <- NA }
        tmp.jname <- custom_grep(all_journals[[i]], "Title", format = "char")
        if (is.null(tmp.jname)) { tmp.jname <- NA }
        # Author specific stuff, requires an extra loop
        tmp.authlist <- custom_grep(all_authors[[i]], "Author", format = "list")
        if ((!is.null(tmp.authlist)) & length(tmp.authlist) > 0) {
          out <- lapply(tmp.authlist, (function(auth){
            tmp.last <- custom_grep(auth, "LastName", format = "char")
            tmp.last <- ifelse(length(tmp.last) > 0, tmp.last, NA)
            #
            tmp.first <- custom_grep(auth, "ForeName", format = "char")
            tmp.first <- ifelse(length(tmp.first) > 0, tmp.first, NA)
            #
            tmp.affi <- custom_grep(auth, "Affiliation", format = "char")
            if (!is.null(tmp.affi)) {
              tmp.affi <- trim_address(tmp.affi[1])
            } else {
              tmp.affi <- NA
            }
            tmp.email <- regexpr("([[:alnum:]]|\\.|\\-\\_){3,200}@([[:alnum:]]|\\.|\\-\\_){3,200}(\\.)([[:alnum:]]){2,6}", 
                                 auth)
            if (tmp.email > 0) {
              tmp.email <- substr(auth, tmp.email[1], tmp.email[1] -1 +
                                attributes(tmp.email)$match.length[1])
            } else {
              tmp.email <- NA
            }
            c(title = tmp.title[1],
              year = tmp.year[1],
              journal = tmp.jname[1],
              firstName = tmp.first[1],
              lastName = tmp.last[1],
              address = tmp.affi[1],
              email = tmp.email[1])
          }))
          out <- suppressWarnings(data.frame(do.call(rbind, out), stringsAsFactors = FALSE))
          #
          # impute addresses
          out <- out[!is.na(out$lastName),]
          ADDKEEP <- !is.na(out$address)
          if(sum(ADDKEEP) > 0 ) {
            ADDKEEP <- which(ADDKEEP)
            RESIDX <- rep(ADDKEEP[1], nrow(out))
            for (zi in 1:nrow(out)) {
              if(zi %in% ADDKEEP)
                RESIDX[zi:nrow(out)] <- zi  
            }
            out$address <- out$address[RESIDX]
          }
          out
          #
        } else {
          NULL
        }
      }))
      xtracted.data <- suppressWarnings(do.call(rbind, data.out))
      message(" Filtering... ", appendLF = FALSE)
      xtracted.data <- xtracted.data[!is.na(xtracted.data$address),
                                     ]
      if (!is.null(affi_regex_exclude)) {
        xtracted.data <- xtracted.data[regexpr(affi_regex_exclude,
                                               xtracted.data$address, ignore.case = TRUE) <
                                         0, ]
      }
      message("done!", appendLF = TRUE)
      xtracted.data
    } else {
      message("An error occurred!")
      NULL
    }
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
