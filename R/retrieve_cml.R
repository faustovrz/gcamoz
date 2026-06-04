## Resolve CML IDs for the 21 CIMMYT maize accessions from Genesys.
## Uses the saved browser/Google token (see R/genesys_auth.R).
## Filtering: the API's `accessionNumber` filter key 500s and most other keys
## are silently ignored; only `_text` (free-text search) actually narrows, so
## we search per accession and exact-match locally on normalized ACCENUMB.

suppressMessages({
  source("R/genesys_auth.R")
  library(tidyverse)
})
genesys_use_saved_token()

target_accessions <- c(
  "CIMMYTMA 29755", "CIMMYTMA 29754", "CIMMYTMA 29753", "CIMMYTMA 29752",
  "CIMMYTMA 29751", "CIMMYTMA 29750", "CIMMYTMA 28038", "CIMMYTMA 28037",
  "CIMMYTMA 28035", "CIMMYTMA 28034", "CIMMYTMA 28033", "CIMMYTMA 27991",
  "CIMMYTMA 27990", "CIMMYTMA 27989", "CIMMYTMA 27988", "CIMMYTMA 27987",
  "CIMMYTMA 27986", "CIMMYTMA 27985", "CIMMYTMA 27984", "CIMMYTMA 27983",
  "CIMMYTMA 27982"
)
normalize_accenumb <- function(x) str_squish(toupper(x))

select <- genesysr:::.fieldsToSelect(
  list("INSTCODE", "ACCENUMB", "DOI", "GENUS", "SPECIES", "ACCENAME")
)
base <- list(
  institute = list(code = list("MEX002")),
  taxonomy  = list(genus = list("Zea"))
)

rows <- map_dfr(target_accessions, function(a) {
  f <- c(base, list("_text" = a))
  d <- tryCatch(
    genesysr:::.list_accessions_page(f, 0, 50, select),
    error = function(e) data.frame()
  )
  if (nrow(d)) d else data.frame()
}) %>% distinct()

cml_pattern <- "CML[[:space:]]?[0-9]+"
rows <- rows %>%
  mutate(
    accenumb_norm = normalize_accenumb(ACCENUMB),
    cml_id = str_remove(str_extract(ACCENAME, cml_pattern), "\\s")
  )

result <- tibble(
  accession_number = target_accessions,
  accenumb_norm    = normalize_accenumb(target_accessions)
) %>%
  left_join(
    rows %>% select(accenumb_norm, cml_id, ACCENAME, DOI),
    by = "accenumb_norm"
  ) %>%
  select(accession_number, cml_id, accession_name = ACCENAME, DOI)

dir.create("output", showWarnings = FALSE)
write_csv(result, "output/cml_ids.csv")
saveRDS(rows, "output/genesys_passport_raw.rds")

cat("\nTotal rows fetched:", nrow(rows),
    "| Matched:", sum(!is.na(result$cml_id)),
    "| Unmatched:", sum(is.na(result$cml_id)), "\n\n")
print(as.data.frame(result), right = FALSE)
