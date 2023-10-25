#' Removes subgenus from a taxon name
#' @description This function removes the subgenus "(subgenus)" from a taxon name.
#' @param df Input dataframe
#' @param name Column containing taxon names
#' @param rank Column containing information on whether a taxon is a subgenus. This column needs to contain the character string "subg" (ignoring case) for any taxa that are subgenera.
#'
#' @export
taxon_subgenus_rm <- function(df, name, rank) {

  name1 <- name2 <- rank1 <- NULL

  colnames(df)[colnames(df) == name] <- "name1"
  colnames(df)[colnames(df) == rank] <- "rank1"

  name_df <- as.data.frame(df[,c("name1", "rank1")])
  colnames(name_df) <- c("name1", "rank1")

  name_df <- name_df %>%
    dplyr::distinct(name1, rank1) %>%
    dplyr::mutate(name2 = dplyr::if_else(!grepl("subg", tolower(rank1)), trimws(gsub("\\([[:alpha:]]*\\)[[:space:]]", "", name1)), trimws(gsub("\\([[:alpha:]]*\\)$", "", name1))))

  name_df[,name] <- name_df[,"name2"]
  name_df[,rank] <- name_df[,"rank1"]

  df <- dplyr::left_join(df, name_df, by = c("name1", "rank1")) %>%
    dplyr::select(-c(name1, name2, rank1))

}
