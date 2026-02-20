# ============================================================
# IsoFrag Custom Functions - Simplified for Phytol Analysis
# ============================================================
# Author: Merve Tomečková Öztoprak
# Description: Core functions needed for MS1 and MS2 phytol isotopomer analysis

# ============= Internal utility functions =============

# Internal function to wrap expressions with error handling
try_catch_all <- function(expr, error, warn = error, newline = TRUE) {
  tryCatch(
    {{ expr }},
    error = function(p) {
      if(newline) cat("\n")
      if (is.function(error)) error(p)
      else stop(error, call. = FALSE)
    },
    warning = function(p) {
      if(newline) cat("\n")
      if (is.function(warn)) warn(p)
      else stop(warn, call. = FALSE)
    }
  )
}

# Group one dataset the same as another (restore original groupings)
group_by_same_groups <- function(target_dataset, source_dataset) {
  target_dataset |> dplyr::group_by(!!!dplyr::groups(source_dataset))
}

# ============= Core flagging functions =============

# Flag satellite peaks based on mass defect
orbi_flag_mass.def_peaks <- function(dataset) {
  
  # safety checks
  cols <- c("filename", "compound", "scan.no", "isotopocule", "ions.incremental", "m_z", "mzMeasured")
  stopifnot(
    "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
    "`dataset` requires required columns" = all(cols %in% names(dataset))
  )
  
  # calculation
  dataset <- try_catch_all(
    dataset |>
      dplyr::group_by(.data$filename, .data$compound, .data$scan.no, .data$isotopocule) |>
      dplyr::mutate(
        massdef.mz = abs(.data$m_z - .data$mzMeasured),
        is_duplicate = n() > 1,
        is_main_mz = .data$ions.incremental == max(.data$ions.incremental)
      ) |>
      dplyr::ungroup(),
    "Error flagging satellite peaks: "
  )
  
  return(dataset)
}

# Flag outliers based on TIC*IT product
orbi_flag_outliers_new <- function(dataset, agc_sd_cutoff = 3) {
  
  # safety checks
  cols <- c("filename", "scan.no", "tic", "it.ms")
  stopifnot(
    "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
    "`dataset` requires columns `filename`, `scan.no`, `tic` and `it.ms`" = all(cols %in% names(dataset))
  )
  
  # calculation
  dataset <- try_catch_all(
    dataset |>
      dplyr::group_by(.data$filename) |>
      dplyr::mutate(
        TICxIT = .data$tic * .data$it.ms,
        med = median(.data$TICxIT, na.rm = TRUE),
        sdv = sd(.data$TICxIT, na.rm = TRUE),
        is_outlier = .data$TICxIT < (.data$med - agc_sd_cutoff * .data$sdv) |
          .data$TICxIT > (.data$med + agc_sd_cutoff * .data$sdv)
      ) |>
      dplyr::ungroup(),
    "Error flagging outliers: "
  )
  
  return(dataset)
}

# Flag weak isotopocules based on scan coverage
orbi_flag_weak_isotopocules <- function(dataset, min_percent = 90) {
  
  # safety checks
  cols <- c("filename", "compound", "scan.no", "isotopocule")
  stopifnot(
    "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
    "`dataset` requires required columns" = all(cols %in% names(dataset))
  )
  
  # calculation
  dataset <- try_catch_all(
    dataset |>
      dplyr::group_by(.data$filename, .data$compound, .data$isotopocule) |>
      dplyr::mutate(
        prcnt.coverage = length(unique(.data$scan.no)) / max(.data$scan.no) * 100,
        is_weak_isotopocule = .data$prcnt.coverage < min_percent
      ) |>
      dplyr::ungroup(),
    "Error flagging weak isotopocules: "
  )
  
  return(dataset)
}

# ============= Base peak and Nio functions =============

# Filter scans missing base peak
orbi_basepeak_filter <- function(clean.data, basepeak_def) {
  
  # safety checks
  stopifnot(
    "need a data frame" = !missing(clean.data) && is.data.frame(clean.data),
    "`basepeak_def` needs to be a single text value" = !missing(basepeak_def) && is.character(basepeak_def) && length(basepeak_def) == 1
  )
  
  # ensure factors
  if("isotopocule" %in% names(clean.data) && !is.factor(clean.data$isotopocule)) {
    clean.data$isotopocule <- factor(clean.data$isotopocule)
  }
  
  # identify `basepeak` for each scan
  df.out <- clean.data |>
    dplyr::group_by(.data$filename, .data$compound, .data$scan.no) |>
    dplyr::filter(any(.data$isotopocule == basepeak_def)) |>
    dplyr::ungroup()
  
  return(df.out)
}

# ============= Summary helper functions =============

# Calculate standard error
calculate_ratios_sem <- function(ratios) {
  stopifnot(
    "no input vector for `ratios` supplied" = !missing(ratios),
    "`ratios` need to be a numeric vector" = is.vector(ratios) && is.numeric(ratios),
    "length of `ratios` needs to be > 1" = length(ratios) > 1L
  )
  
  stats::sd(ratios, na.rm = TRUE) / sqrt(length(ratios))
}

# ============= Fractional abundance function =============

# Calculate fractional abundances
f.abundance <- function(x) {
  data <- x %>%
    dplyr::mutate(element = stringr::str_extract(.data$isotopocule, "(?<=\\d)\\p{L}+")) %>%
    dplyr::group_by(.data$filename, .data$compound, .data$isotopocule) %>%
    dplyr::mutate(n.total.scans = dplyr::n()) %>%
    dplyr::ungroup()
  
  data_out <- data[order(data$filename, data$compound, data$isotopocule, data$scan.no), ]
  
  data_out <- data_out %>%
    # Calculate relative abundance of 2H
    dplyr::group_by(.data$filename, .data$compound) %>%
    dplyr::mutate(rel.ab.2H = length(which(.data$isotopocule == "2H")) / length(unique(.data$scan.no)) * 100) %>%
    dplyr::ungroup() %>%
    # Calculate Nio
    dplyr::mutate(Nio = ((.data$intensity / .data$peakNoise) * (4.4 / 1)) * 
                    ((120000 / .data$resolution)^0.5) * (.data$microscans^0.5)) %>%
    # Complete missing isotopocules
    tidyr::complete(nesting(.data$filename, .data$compound, .data$scan.no),
                    .data$isotopocule,
                    fill = list(Nio = 0),
                    explicit = FALSE) %>%
    dplyr::mutate(
      simulated = is.na(.data$intensity),
      element = dplyr::case_when(
        .data$simulated == TRUE & .data$isotopocule == "0U" ~ "U",
        .data$simulated == TRUE & .data$isotopocule == "13C" ~ "C",
        .data$simulated == TRUE & .data$isotopocule == "2H" ~ "H",
        .data$simulated == TRUE & .data$isotopocule == "17O" ~ "O",
        TRUE ~ .data$element
      ),
      Nio = ifelse(.data$simulated == TRUE & .data$element == "U", 1, .data$Nio)
    ) %>%
    # Extract n.H and n.C from compound formula
    dplyr::mutate(
      n.H = dplyr::case_when(
        .data$compound %in% c("C5H7") ~ 7,
        .data$compound %in% c("C5H9", "C5H9O", "C6H9") ~ 9,
        .data$compound %in% c("C5H11", "C6H11", "C6H11O", "C7H11") ~ 11,
        .data$compound %in% c("C6H13", "C7H13", "C7H13O", "C8H13") ~ 13,
        .data$compound %in% c("C7H15", "C8H15", "C8H15O", "C9H15") ~ 15,
        .data$compound %in% c("C8H17", "C9H17", "C9H17O", "C10H17") ~ 17,
        .data$compound %in% c("C9H19", "C10H19", "C11H19") ~ 19,
        .data$compound %in% c("C11H21", "C12H21") ~ 21,
        .data$compound == "C12H23" ~ 23
      ),
      n.C = dplyr::case_when(
        .data$compound %in% c("C5H7", "C5H9", "C5H9O", "C5H11") ~ 5,
        .data$compound %in% c("C6H9", "C6H11", "C6H13", "C6H11O") ~ 6,
        .data$compound %in% c("C7H11", "C7H13", "C7H15", "C7H13O") ~ 7,
        .data$compound %in% c("C8H17", "C8H13", "C8H15", "C8H15O") ~ 8,
        .data$compound %in% c("C9H15", "C9H17", "C9H19", "C9H17O") ~ 9,
        .data$compound %in% c("C10H17", "C10H19") ~ 10,
        .data$compound %in% c("C11H19", "C11H21") ~ 11,
        .data$compound %in% c("C12H21", "C12H23") ~ 12
      )
    ) %>%
    # Calculate fractional abundances
    dplyr::group_by(.data$filename, .data$compound, .data$scan.no) %>%
    dplyr::mutate(
      fa_denom_sum_Nio_noH = ifelse(.data$isotopocule != "2H", sum(.data$Nio), NA),
      fa_Nio_noH = .data$Nio / .data$fa_denom_sum_Nio_noH,
      fa_denom_sum_Nio = sum(.data$Nio),
      fa_Nio = .data$Nio / .data$fa_denom_sum_Nio
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$filename, .data$compound, .data$isotopocule, .data$element, 
                  .data$n.H, .data$n.C, .data$scan.no, .data$time.min, .data$m_z,
                  .data$n.total.scans, .data$tic, .data$it.ms, .data$intensity,
                  .data$ions.incremental, .data$Nio, .data$fa_Nio, .data$fa_Nio_noH,
                  .data$rel.ab.2H, .data$simulated) %>%
    dplyr::arrange(.data$filename, .data$isotopocule, .data$n.C, .data$n.H)
  
  return(data_out)
}