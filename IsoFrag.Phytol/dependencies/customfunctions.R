#custom functions----
#author: Merve Tomečková Öztoprak

#base functions are copied from isoorbi package
# Internal utility functions =============

#@title Internal function to calculate standard error
# @description The function `calculate_ratios_sem()` computes a regular standard error.
# @keywords internal
# @param ratios A numeric vector used to calculate a standard error
# @return The calculated standard error
calculate_ratios_sem <- function(ratios) {
  
  # safety checks
  stopifnot(
    "no input vector for `ratios` supplied" = !missing(ratios),
    "`ratios` need to be provided as a numeric vector" = is.vector(ratios) && is.numeric(ratios),
    "length of `ratios` needs to be > 1, cannot calculate SEM with a single value" = length(ratios) > 1L
  )
  
  # calculation
  try_catch_all(
    stats::sd(ratios) / sqrt(length(ratios)),
    "something went wrong calculating the standard error: "
  )
}

# @title Internal function to calculate geometric mean
# @description  The function `calculate_ratios_gmean()` is used to calculate geometric means.
# @keywords internal
# @param ratios A numeric vector of ratios used to calculate the geometric mean
# @return The calculated geometric mean
calculate_ratios_gmean <- function(ratios) {
  
  # safety checks
  stopifnot(
    "no input vector for `ratios` supplied" = !missing(ratios),
    "`ratios` need to be provided as a numeric vector" = is.vector(ratios) && is.numeric(ratios)
  )
  
  # calculation
  try_catch_all(
    exp(mean(log(ratios))),
    "something went wrong calculating the geometic mean: "
  )
}

# internal function to ensure dataset columns are factors (if they exist)
# will warn the user about the transformation
factorize_dataset <- function(dataset, cols = c()) {
  for (col in cols) {
    if (col %in% names(dataset) && !is.factor(dataset[[col]])) {
      sprintf("column `%s` was turned into a factor", col) |> message()
      dataset[[col]] <- factor_in_order(dataset[[col]])
    }
  }
  return(dataset)
}

# group if exists
group_if_exists <- function(dataset, cols, add = TRUE) {
  for (col in cols) {
    if (col %in% names(dataset))
      dataset <- dataset |> group_by(!!sym(col), .add = add)
  }
  return(dataset)
}

# group one dataset the same as another (use to restore original groupings)
group_by_same_groups <- function(target_dataset, source_dataset) {
  target_dataset |> dplyr::group_by(!!!dplyr::groups(source_dataset))
}

# count distinct column respecting existing grouping
count_grouped_distinct <- function(dataset, column) {
  dataset |> dplyr::select(!!!c(dplyr::groups(dataset), column)) |>
    dplyr::distinct() |> dplyr::count(!!sym(column)) |>
    dplyr::pull(.data$n) |> sum()
}

# internal function to wrap around expressions that should throw errors whenever anything unexpected happens
# error/warn can be text or function
try_catch_all <- function(expr, error, warn = error, newline = TRUE) {
  tryCatch(
    {{ expr }},
    error = function(p) {
      if(newline) cat("\n")
      if (is_function(error)) error(p)
      else abort(error, parent = p)
    },
    warning = function(p) {
      if(newline) cat("\n")
      if (is_function(warn)) warn(p)
      else abort(warn, parent = p)
    }
  )
}

# print out info start message
message_start <- function(...) {
  message_wrap(..., exdent = 3, appendLF = !interactive())
  return(Sys.time())
}

# print out info end message
message_finish <- function(..., start_time = NULL, pre = if (!interactive()) "..." else "",  indent = if (interactive()) 0 else 3) {
  end_time <- Sys.time()
  time_info <-
    if (!is.null(start_time)) sprintf(" in %.2f seconds.", end_time - start_time)
  else ""
  message_wrap(pre, ..., time_info, indent = indent, exdent = 3)
}

# print out standalone info message
message_standalone <- function(..., start_time = NULL) {
  message_finish(..., start_time = start_time, pre = "", indent = 0)
}

# print out a message that wraps in non-interactive mode (e.g. in notebooks)
message_wrap <- function (..., appendLF = TRUE, width = if (!interactive()) options()$width - 2 else NA, indent = 0, exdent = 0) {
  if (is.na(width)) {
    message(..., appendLF = appendLF)
  } else {
    original <- paste0(..., collapse = "")
    original |>
      strwrap(width = width, indent = indent, exdent = exdent) |>
      paste(collapse = "\n") |>
      message(appendLF = appendLF)
  }
  invisible()
}

# internal function to filter for specific isotopocules
filter_isotopocules <- function(dataset, isotopocules, allow_all = TRUE) {
  dataset <- dataset |> factorize_dataset("isotopocules")
  if (allow_all && length(isotopocules) == 0L)
    isotopocules <- levels(dataset$isotopocule)
  missing_isotopocules <- !isotopocules %in% levels(dataset$isotopocule)
  if (sum(missing_isotopocules) > 0L) {
    sprintf("not all `isotopocules` are in the dataset, missing '%s'. Available: '%s'",
            paste(isotopocules[missing_isotopocules], collapse = "', '"),
            paste(levels(dataset$isotopocule), collapse = "', '")) |>
      warn()
  }
  isotopocules <- isotopocules[!missing_isotopocules]
  if (length(isotopocules) == 0L)
    abort("none of the provided `isotopocules` are in the dataset")
  
  # plot dataset
  dataset |>
    dplyr::filter(.data$isotopocule %in% isotopocules) |>
    droplevels()
}

# internal function for the facet_wrap
# decides whether to wrap by filename, compound or both filename and compound
# depending on if either has more than 1 value
dynamic_wrap <- function(plot, scales = "free_x") {
  dataset <- plot$data |> factorize_dataset(c("filename", "compound"))
  n_files <- length(levels(dataset$filename))
  n_compounds <- length(levels(dataset$compound))
  if (n_files > 1L && n_compounds > 1L) {
    plot <- plot + ggplot2::facet_wrap(~.data$filename + .data$compound, scales = scales)
  } else if (n_compounds > 1L) {
    plot <- plot + ggplot2::facet_wrap(~.data$compound, scales = scales)
  } else if (n_files > 1L) {
    plot <- plot + ggplot2::facet_wrap(~.data$filename, scales = scales)
  }
  return(plot)
}

# y axis as log, pseudo-log, or continuous
dynamic_y_scale <- function(plot,  y_scale = c("raw", "linear", "pseudo-log", "log"), sci_labels = FALSE, breaks = scales::pretty_breaks(5)) {
  y_scale <- arg_match(y_scale)
  labeler <- if(sci_labels) label_scientific_log() else identity
  if (y_scale == "log") {
    plot <- plot + 
      ggplot2::scale_y_log10(label = labeler) +
      ggplot2::annotation_logticks(sides = "l")
  } else if (y_scale == "pseudo-log") {
    plot <- plot +
      ggplot2::scale_y_continuous(
        trans = scales::pseudo_log_trans(),
        breaks = breaks,
        labels = labeler
      )
  } else if (y_scale == "linear"){
    plot <- plot +
      ggplot2::scale_y_continuous(labels = labeler)
  }
  return(plot)
}

# internal function to nicely format log scales
#' @importFrom stats na.omit
label_scientific_log <- function() {
  parser1 <- scales::label_scientific()
  parser2 <- scales::label_parse()
  parser3 <- scales::label_log()
  function(x) {
    needs_decimal <- any((log10(na.omit(x)) %% 1) > 0)
    if (needs_decimal) {
      parsed_x <- x |>
        parser1()
      out <- sub("e\\+?", " %.% 10^", parsed_x)
      out <- out |>
        parser2()
    } else {
      out <- parser3(x)
    }
    out[x == 0.0] <- 0
    return(out)
  }
}

#' Default isoorbi plotting theme
#' @return ggplot theme object
#' @param text_size a font size for text
#' @param facet_text_size a font size for facet text
#' @export
orbi_default_theme <- function(text_size = 16, facet_text_size = 20) {
  ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = text_size),
      strip.text = ggplot2::element_text(size = facet_text_size),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    )
}

# =========== my functions ==============

# import function
import_data <- function(filepath) {
  raw.files <- tibble(filepath = list.files(filepath, pattern = "*.isox", full.names = TRUE))
  
  raw.files %>%
    mutate(filename = trimws(filepath, whitespace = ".*\\/")) %>%
    relocate(filename) %>%
    rowwise() %>%
    do(., orbi_read_isox(file = .$filepath))
}

#quick summary on imported data
process_raw_data <- function(raw.data) {
  raw.data %>%
    rename(isotopolog = isotopocule) %>%
    group_by(filename) %>%
    summarise(
      no.scans = length(unique(scan.no)),
      time = max(time.min),
      compound = length(unique(compound)),
      isotopolog = length(unique(isotopolog)),
      m.resolution = mean(resolution),
      m.microscans = mean(microscans),
      m.AGC = mean(agcTarget)
    )
}

# Data Filtering and Quality Control
filter_data <- function(data) {
  # Filter isotopocules and compounds
  filtered <- data %>%
    filter(isotopocule %in% c("M0", "13C", "2H")) %>%
    filter(!compound %in% c("C20H37", "C20H39O", "C20H40O"))
  
  # Flag and remove outliers
  flag_data <- filtered %>%
    orbi_flag_outliers_new(agc_sd_cutoff = 2) %>%
    orbi_flag_mass.def_peaks() %>%
    orbi_flag_weak_isotopocules(min_percent = 90)
}

# Isotope Abundance Calculations 
calculate_abundances <- function(dataset) {
  dataset %>%
    mutate(
      element = str_extract(isotopocule, "(?<=\\d)\\p{L}+"),
      Nio = ((intensity/peakNoise)*(4.4/1))*((120000/resolution)^0.5)*(microscans^0.5),
      n.C = as.numeric(str_extract(compound, "(?<=C)\\d+")),
      n.H = as.numeric(str_extract(compound, "(?<=H)\\d+"))
    ) %>%
    select(filename, compound, n.C, n.H, isotopocule, scan.no, Nio, is_weak_isotopocule)
}

# Naive approach fractional abundance calculation
naive_approach <- function(data) {
  data %>%
    complete(nesting(filename, compound, scan.no), isotopocule, 
             fill = list(Nio = NA), explicit = FALSE) %>%
    group_by(filename, compound, isotopocule) %>%
    split(.$isotopocule) %>%
    as.data.frame() %>%
    rename(filename = M0.filename, compound = M0.compound, scan.no = M0.scan.no) %>%
    select(filename, compound, scan.no, all_of(grep("Nio$", names(.), value = TRUE))) %>%
    mutate(ro_denom.sum = rowSums(select(., ends_with(".Nio")), na.rm = TRUE)) %>%
    mutate(across(all_of(grep("(13C|2H|M0).*Nio", names(.), value = TRUE)), 
                  ~ ./ro_denom.sum, 
                  .names = "ro.{gsub('X', '', sub('.Nio$', '', col))}")) %>%
    complete(filename, compound) %>%
    mutate(sum.ro = rowSums(across(starts_with("ro.")), na.rm = TRUE)) %>%
    group_by(filename, compound) %>%
    mutate(shotnoise.H = (((1/sum(na.omit(X2H.Nio)) + (1/sum(na.omit(M0.Nio)))))^0.5))
}

# Naive summary calculation
calculate_stats <- function(data, method_name, iso_pattern = "ro.", elements = c("M0", "13C", "2H")) {
  # Validate inputs
  stopifnot(is.data.frame(data))
  stopifnot(is.character(method_name), length(method_name) == 1)
  required_cols <- c("filename", "compound", "scan.no", "shotnoise.H", "ro.13C")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  # Process data
  result <- data %>%
    group_by(filename, compound,shotnoise.H) %>%
    summarise(
      method = method_name,
      across(starts_with(iso_pattern), 
             list(
               m = ~ mean(.x, na.rm = TRUE),
               ae = ~ (sd(.x, na.rm = TRUE) / sqrt(n()))
             ),
             .names = "{col}.{.fn}"
      ),
      .groups = "drop"
    ) %>%
    mutate(
      across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)),
      sum.ro = rowSums(across(ends_with(".m")), na.rm = TRUE),
      across(ends_with(".ae"),
             ~ ifelse(. <= 2 * shotnoise.H, "pass", "fail"),
             .names = "{col}.pass"),
      n.C = as.numeric(str_extract(compound, "(?<=C)\\d+")),
      n.H = as.numeric(str_extract(compound, "(?<=H)\\d+"))
    ) %>%
    filter(!sum.ro == 0)
  # Select element-specific columns dynamically
  element_cols <- map(elements, ~ paste0(iso_pattern, .x, c(".m", ".ae"))) %>% 
    flatten_chr()
  result %>%
    select(filename, method, compound, n.H, n.C, 
           all_of(element_cols), sum.ro, shotnoise.H,
           ro.M0.ae, ro.13C.ae, ro.2H.ae, ro.13C.ae.pass)
}

# replicate measurement summary per ID ----
summary_ID <- function(s_ro.corr){
  # Process data
  data = s_ro.corr%>%
    group_by(ID,compound,n.C,n.H)%>%
    summarize(
      across(starts_with("ro.")& !ends_with("ae"), 
             list(
               m = ~ mean(.x, na.rm = TRUE),  # Calculate mean
               #relative standard error in replicate measurements! 
               rse = ~ (sd(.x, na.rm = TRUE) / sqrt(n())) *  (1/mean(.x, na.rm = TRUE)),
               error= ~sd(.x, na.rm=TRUE) / sqrt(length(.x))
             ),
             .names = "{.fn}.{.col}"),
      across(ends_with("ae"), 
             list(
               m.ae = ~ mean(.x, na.rm = TRUE),  # Calculate mean
               #error here is the error in replicate measurements! 
               prop.ae = ~ sqrt(sum(.x^2)/ length(.x)^2)
             ),
             .names = "{.fn}.{.col}"
      )
    )%>%
    mutate(DoU=((2*n.C+2)-(n.H+1))/2)
}

# Interitance matrix sum model
process_sum_model <- function(df_group, fragment_definitions) {
  observed_values <- df_group$ro.13C.corr
  names(observed_values) <- df_group$compound
  
  valid_fragments <- intersect(names(fragment_definitions), df_group$compound)
  fragment_definitions <- fragment_definitions[valid_fragments]
  observed_values <- observed_values[valid_fragments]
  
  A_sum <- matrix(0, nrow = length(valid_fragments), ncol = 20,
                  dimnames = list(valid_fragments, carbon_positions))
  for (fragment in valid_fragments) {
    A_sum[fragment, fragment_definitions[[fragment]]] <- 1
  }
  
  solution_sum <- lm(observed_values ~ A_sum + 0)
  position_specific_13C_sum <- coef(solution_sum)
  names(position_specific_13C_sum) <- carbon_positions
  
  position_specific_13C_sum_noNA <- position_specific_13C_sum
  position_specific_13C_sum_noNA[is.na(position_specific_13C_sum_noNA)] <- 0
  
  # Post hoc constraints
  position_specific_13C_sum_noNA[1:3] <- position_specific_13C_sum[1] / 4
  position_specific_13C_sum_noNA[17] <- position_specific_13C_sum[1] / 4
  position_specific_13C_sum_noNA[14:16] <- position_specific_13C_sum[14] / 4
  position_specific_13C_sum_noNA[20] <- position_specific_13C_sum[14] / 4
  
  return(position_specific_13C_sum_noNA)
}

# Interitance matrix sum model standard sample
process_sum_model_STD <- function(df_group, fragment_definitions) {
  observed_values <- df_group$delta_calc
  names(observed_values) <- df_group$compound
  
  valid_fragments <- intersect(names(fragment_definitions), df_group$compound)
  fragment_definitions <- fragment_definitions[valid_fragments]
  observed_values <- observed_values[valid_fragments]
  
  A_sum <- matrix(0, nrow = length(valid_fragments), ncol = 20,
                  dimnames = list(valid_fragments, carbon_positions))
  for (fragment in valid_fragments) {
    A_sum[fragment, fragment_definitions[[fragment]]] <- 1
  }
  
  solution_sum <- lm(observed_values ~ A_sum + 0)
  position_specific_13C_sum <- coef(solution_sum)
  names(position_specific_13C_sum) <- carbon_positions
  
  position_specific_13C_sum_noNA <- position_specific_13C_sum
  position_specific_13C_sum_noNA[is.na(position_specific_13C_sum_noNA)] <- 0
  
  # Post hoc constraints
  position_specific_13C_sum_noNA[1:3] <- position_specific_13C_sum[1] / 4
  position_specific_13C_sum_noNA[17] <- position_specific_13C_sum[1] / 4
  position_specific_13C_sum_noNA[14:16] <- position_specific_13C_sum[14] / 4
  position_specific_13C_sum_noNA[20] <- position_specific_13C_sum[14] / 4
  
  return(position_specific_13C_sum_noNA)
}


# Define which carbon positions belong to each fragment Model
#-------
fragment_definitions <- list(
  C4H7 = c("C16", "C15", "C14", "C20"),
  C4H9 = c("C16", "C15", "C14", "C20"),
  C4H7O = c("C1", "C2", "C3", "C17"),
  C5H7 = c("C1", "C2", "C3", "C17", "C4"),
  C5H9 = c("C16", "C15", "C14", "C20","C13"),
  C5H11 = c("C16", "C15", "C14", "C20","C13"),
  C5H9O = c("C1", "C2", "C3", "C17", "C4"),
  C6H9 = c("C1", "C2", "C3", "C17", "C4", "C5"),
  C6H11 = c("C16", "C15", "C14", "C20", "C13", "C12"),
  C6H13 = c("C16", "C15", "C14", "C20", "C13", "C12"),
  C6H11O = c("C1", "C2", "C3", "C17", "C4", "C5"),
  C7H11 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6"),
  C7H13 = c("C16", "C15", "C14", "C13", "C12", "C11", "C20"),
  C7H15 = c("C16", "C15", "C14", "C13", "C12", "C11", "C20"),
  C7H13O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6"),
  C8H13 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7"),
  C8H15 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20"),
  C8H17 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20"),
  C8H15O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7"),
  C9H15 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18"),
  C9H17 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10"),
  C9H19 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10"),
  C9H17O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18"),
  C10H17 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8"),
  C10H19 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9"),
  C10H21 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9"),
  C10H19O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8"),
  C11H19 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9"),
  C11H21 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8"),
  C11H23 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8"),
  C11H21O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9"),
  C12H21 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9", "C10"),
  C12H23 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8", "C7"),
  C12H25 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8", "C7"),
  C12H23O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9", "C10"),
  C13H23 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9", "C10", "C11"),
  C13H25 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8", "C7", "C18"),
  C13H27 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8", "C7", "C18"),
  C13H25O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9", "C10", "C11")
)
#-----

# alternative sattelite peak identification based on minimum relative mass defect difference 
orbi_flag_mass.def_peaks <- function(dataset) {
   
   # safety checks
   cols <- c("filepath", "filename", "compound", "scan.no", "time.min", "isotopocule", "ions.incremental", "tic", "it.ms","m_z")
   stopifnot(
     "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
     "`dataset` requires columns `filepath`, `filename`, `compound`, `scan.no`, `time.min`, `isotopocule`, `ions.incremental`, `tic`, `it.ms` and 'm_z'" =
       all(cols %in% names(dataset))
   )
   
   # info
   start_time <- message_start("orbi_flag_mass.def_peaks() is flagging minor signals (identified by reltive mass defect)... ")
   
   # calculation
   dataset <-
     try_catch_all(
       dataset |>
         dplyr::group_by(
           .data$filename,
           .data$compound,
           .data$scan.no,
           .data$isotopocule
         ) |>
         dplyr::mutate(
           massdef.mz = abs(.data$m_z-.data$mzMeasured),
           is_duplicate = n()>1,
           is_main_mz=.data$ions.incremental == max(.data$ions.incremental),
           ) |># calculate mean mass measured for non duplicate scans for each compound
         dplyr::ungroup(),
       "something went wrong tying to flag satellite peaks: "
     )
   
   # info
   sat_mass.def_peaks<- sum(dataset$is_main_mz == FALSE & dataset$is_duplicate == TRUE)
   
   isotopocules <- dataset |> 
     dplyr::summarise(has_mass.deff.sat_peaks = any(.data$is_main_mz == FALSE) & any(dataset$is_duplicate == TRUE), .by = "isotopocule") |>
     dplyr::filter(.data$has_mass.deff.sat_peaks) |>
     dplyr::pull(.data$isotopocule)
   
   if(sat_mass.def_peaks > 0){
     sprintf(
       "flagged %d/%d peaks in %d isotopocules (\"%s\") as satellite peaks (%.1f%%) identified by relative mass defect",
       sat_mass.def_peaks, nrow(dataset), length(isotopocules), 
       paste(isotopocules, collapse = "\", \""),
       100 * sat_mass.def_peaks/nrow(dataset)) |>
       message_finish(start_time = start_time)
   } else {
     message_finish("confirmed there are no satellite peaks", start_time = start_time)
   }
   return(dataset)
 }
orbi_flag_mass.def_peaks2 <- function(dataset) {
  
  # safety checks
  cols <- c("filepath", "filename", "compound", "scan.no", "time.min", "isotopocule", "ions.incremental", "tic", "it.ms","m_z")
  stopifnot(
    "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
    "`dataset` requires columns `filepath`, `filename`, `compound`, `scan.no`, `time.min`, `isotopocule`, `ions.incremental`, `tic`, `it.ms` and 'm_z'" =
      all(cols %in% names(dataset))
  )
  
  # info
  start_time <- message_start("orbi_flag_mass.def_peaks() is flagging minor signals (identified by reltive mass defect)... ")
  
  # calculation
  dataset <-
    try_catch_all(
      dataset |>
        dplyr::group_by(
          .data$filename,
          .data$compound,
          .data$scan.no,
          .data$isotopocule
        ) |>
        dplyr::mutate(
          massdef.mz = abs(.data$m_z-.data$mzMeasured),
          is_duplicate = n()>1,
          is_main_mz=.data$massdef.mz == min(.data$massdef.mz),
        ) |># calculate mean mass measured for non duplicate scans for each compound
        dplyr::ungroup(),
      "something went wrong tying to flag satellite peaks: "
    )
  
  # info
  sat_mass.def_peaks<- sum(dataset$is_main_mz == FALSE & dataset$is_duplicate == TRUE)
  
  isotopocules <- dataset |> 
    dplyr::summarise(has_mass.deff.sat_peaks = any(.data$is_main_mz == FALSE) & any(dataset$is_duplicate == TRUE), .by = "isotopocule") |>
    dplyr::filter(.data$has_mass.deff.sat_peaks) |>
    dplyr::pull(.data$isotopocule)
  
  if(sat_mass.def_peaks > 0){
    sprintf(
      "flagged %d/%d peaks in %d isotopocules (\"%s\") as satellite peaks (%.1f%%) identified by relative mass defect",
      sat_mass.def_peaks, nrow(dataset), length(isotopocules), 
      paste(isotopocules, collapse = "\", \""),
      100 * sat_mass.def_peaks/nrow(dataset)) |>
      message_finish(start_time = start_time)
  } else {
    message_finish("confirmed there are no satellite peaks", start_time = start_time)
  }
  return(dataset)
}

# alternative orbi_filter_isox to work with separate csv file in folder to trim acquisition time
alt_orbi_filter_isox <-
  function(dataset,
           filenames = NULL,
           compounds = NULL,
           isotopocules = NULL
           ) {
    
    # safety checks
    cols <- c("filename", "compound", "isotopocule", "time.min", "time_min", "time_max")
    stopifnot(
      "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
      "`dataset` requires columns `filepath`, `filename`, `compound`, `scan.no`, `tic`, `it.ms`, `time_min` and `time_max`" =
        all(cols %in% names(dataset)),
      "`filenames` must be a vector of filenames (or NULL)" = is.null(filenames) || is_character(filenames),
      "`compounds` must be a vector of compounds (or NULL)" = is.null(compounds) || is_character(compounds),
      "`isotopocules` must be a vector of isotopocules (or NULL)" = is.null(isotopocules) || is_character(isotopocules)
    )
    
    # info
    filters <- c()
    if(!is.null(filenames))
      filters <- c(filters, sprintf("filenames (%s)", paste(filenames, collapse = ", ")))
    if(!is.null(compounds))
      filters <- c(filters, sprintf("compounds (%s)", paste(compounds, collapse = ", ")))
    if(!is.null(isotopocules))
      filters <- c(filters, sprintf("isotopocules (%s)", paste(isotopocules, collapse = ", ")))
    
    # info
    n_row_start <- nrow(dataset)
    start_time <-
      sprintf(
        "alt_orbi_filter_isox() is filtering the dataset by %s... ", paste(filters, collapse = ", ")
      ) |>
      message_start()
    
    # filtering
    dataset <-
      try_catch_all({
        # file: filenames
        if (!is.null(filenames))
          dataset <- dataset |> dplyr::filter(.data$filename %in% !!filenames)
        
        # filter: compounds
        if (!is.null(compounds))
          dataset <- dataset |> dplyr::filter(.data$compound %in% !!compounds)
        
        # filter: isotopocules
        if (!is.null(isotopocules))
          dataset <- dataset |> dplyr::filter(.data$isotopocule %in% !!isotopocules)
        
        # filter: time_min
        dataset <- dataset |> dplyr::filter(.data$time.min >= .data$time_min)
        
        # filter: time_max
        dataset <- dataset |> dplyr::filter(.data$time.min <= .data$time_max)
        
        # return
        dataset
      },
      "something went wrong trying to filter dataset: "
      )
    
    # info
    sprintf(
      "removed %d/%d rows (%.1f%%)",
      n_row_start - nrow(dataset), n_row_start, (n_row_start - nrow(dataset))/n_row_start * 100) |>
      message_finish(start_time = start_time)
    
    # return
    return(dataset |> droplevels())
  }
# alternative sattelite peak plotting function based on minimum relative mass defect difference 
orbi_plot_mass.def_peaks <- function(
    dataset, isotopocules = c(), x = c("scan.no", "time.min"), x_breaks = scales::breaks_pretty(5),
    y_scale = c("log", "pseudo-log", "linear", "raw"), y_scale_sci_labels = TRUE,
    colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), 
    color_scale = scale_color_manual(values = colors)) {
  
  # safety checks
  cols <- c("filename", "compound", "scan.no", "time.min", "isotopocule", "ions.incremental")
  stopifnot(
    "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
    "`isotopocules` has to be a character vector if provided" = length(isotopocules) == 0L || is_character(isotopocules),
    "`dataset` requires columns `filename`, `compound`, `scan.no`, `time.min`, `isotopocule`, `ions.incremental`" =
      all(cols %in% names(dataset)),
    "`dataset` requires column `is_notmain_mz` - make sure to run `orbi_flag_mass.def_peaks()` first" = "is_notmain_mz" %in% names(dataset)
  )
  x_column <- arg_match(x)
  y_scale <- arg_match(y_scale)
  
  # prepare dataset
  plot_df <- dataset |>
    factorize_dataset(c("filename", "compound", "isotopocule")) |>
    filter_isotopocules(isotopocules)
  
  # make plot
  plot <- plot_df |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = !!sym(x_column), y = .data$ions.incremental,
      color = .data$isotopocule) +
    ggplot2::geom_line(
      data = function(df) dplyr::filter(df, !.data$is_notmain_mz),
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = function(df) dplyr::filter(df, .data$is_notmain_mz) |>
        dplyr::mutate(flagged = "not min. massdef."),
      map = ggplot2::aes(shape = .data$flagged)
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    ggplot2::scale_shape_manual(values = 17) +
    {{ color_scale }} +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(shape = NA), order = 1)
    ) +
    orbi_default_theme()
  
  # return
  plot |>
    dynamic_y_scale(y_scale, sci_labels = y_scale_sci_labels) |>
    dynamic_wrap()
}

#removes all scans where base peak is missing
orbi_basepeak_filter <- function(clean.data, basepeak_def){
  # safety checks
  stopifnot(
    "need a data frame" = !missing(clean.data) && is.data.frame(clean.data),
    "`clean.data` requires columns `filename`, `compound`, `scan.no`, `isotopocule`, and `ions.incremental`" =
      all(c("filename", "compound", "scan.no", "isotopocule", "ions.incremental") %in% names(clean.data)),
    "`basepeak_def` needs to be a single text value identifying the isotopocule to use as the basepeak" =
      !missing(basepeak_def) && rlang::is_scalar_character(basepeak_def),
    "`basepeak_def` is not an isotopocule in the data frame" = basepeak_def %in% levels(clean.data$isotopocule)
  )
  
  # ensure factors
  clean.data <- clean.data |> factorize_dataset("isotopocule")
  
  # info message
  start_time <-
    sprintf(
      "orbi_basepeak_filter is removing scans of compunds without '%s' isotopocule ... ",
      basepeak_def) |>
    message_start()
  
  # identify `basepeak` for each scan
  df.out <- clean.data |>
    group_by(.data$filename, .data$compound, .data$scan.no) |>
    filter(any(isotopocule == basepeak_def))
  
  # return
  df.out |>
    dplyr::ungroup()
  
  # info message
  sprintf(
    "removed %d (%.1f%%) scans missing base peak '%s' isotopocule ",
    abs(nrow(clean.data)-nrow(df.out)),(abs(nrow(clean.data)-nrow(df.out))/nrow(clean.data))*100,basepeak_def
  ) |> 
    message_finish(start_time = start_time)
  
  # return
  return(df.out)
}

#outlier filter
orbi_flag_outliers_new <- function(dataset, agc_fold_cutoff = NA_real_, agc_window = c(), agc_sd_cutoff = NA_real_) {
  
  # safety checks
  cols <- c("filename", "compound", "scan.no", "tic", "it.ms")
  stopifnot(
    "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
    "`dataset` requires columns `filename`, `compound`, `scan.no`, `tic` and `it.ms`" =
      all(cols %in% names(dataset)),
    "if provided, `agc_fold_cutoff` needs to be a single number" = (length(agc_fold_cutoff) == 1L && is.na(agc_fold_cutoff)) || is_scalar_double(agc_fold_cutoff),
    "if provided, `agc_window` needs to be a vector of two numbers (low and high filter) between 0 and 100" = length(agc_window) == 0L || (is.numeric(agc_window) && length(agc_window) == 2L && agc_window[1] < agc_window[2]),
    "if provided, `agc_sd_cutoff` needs to be a single number" = (length(agc_sd_cutoff) == 1L && is.na(agc_sd_cutoff)) || is_scalar_double(agc_sd_cutoff)
  )
  
  # check filter to apply
  method <- c(
    agc_window = length(agc_window) == 2L,
    agc_fold_cutoff = !is.na(agc_fold_cutoff),
    agc_sd_cutoff = !is.na(agc_sd_cutoff)
  )
  if (sum(method) > 1L) {
    sprintf("can only use one method at a time, please call this function sequentially for each of these parameters: '%s'",
            paste(names(method)[method], collapse = "', '")) |>
      abort()
  } else if (sum(method) == 0L) {
    sprintf("need to define at least one of these parameters for identifying outliers: '%s'",
            paste(names(method), collapse = "', '")) |>
      abort()
  }
  method <- names(method)[method]
  
  # method message
  method_msg <- ""
  method_type <- ""
  if (method == "agc_window") {
    method_msg <- sprintf(
      "%s %% of scans with the lowest and above %s %% of scans with the highest number of ions (`tic` * `it.ms`) in the Orbitrap analyzer",
      agc_window[1], agc_window[2]
    )
    method_type <- sprintf("AGC window (%s to %s %%)", agc_window[1], agc_window[2])
  } else if (method == "agc_fold_cutoff") {
    method_msg <- sprintf(
      "scans below 1/%s and above %s times the average number of ions (`tic` * `it.ms`) in the Orbitrap analyzer",
      agc_fold_cutoff, agc_fold_cutoff
    )
    method_type <- sprintf("%s fold AGC cutoff", agc_fold_cutoff)
  } else if (method == "agc_sd_cutoff") {
    method_msg <- sprintf(
      "scans where `TICxIT` is more than %s times the standard deviation from the median of all scans in a `filename`",
      agc_sd_cutoff
    )
    method_type <- sprintf("%s SD cutoff from median", agc_sd_cutoff)
  }
  
  # optional groupings
  dataset_out <- dataset |> group_if_exists(c("filename", "block", "segment", "injection"), add = TRUE)
  
  # info
  start_time <-
    sprintf(
      "orbi_flag_outliers() is flagging the %s in %d data group(s) (based on '%s')... ",
      method_msg,
      dplyr::n_groups(dataset_out), paste(dplyr::group_vars(dataset_out), collapse = "', '")) |>
    message_start()
  n_scans <- dataset_out |> count_grouped_distinct("scan.no")
  
  # calculation
  dataset_out <- try_catch_all(
    if (method == "agc_window") {
      dataset_out |>
        dplyr::mutate(
          TICxIT = .data$tic * .data$it.ms,
          is_new_outlier =
            .data$TICxIT < stats::quantile(.data$TICxIT, agc_window[1] / 100) |
            .data$TICxIT > stats::quantile(.data$TICxIT, agc_window[2] / 100)
        ) |>
        dplyr::select(-"TICxIT")
    } else if (method == "agc_fold_cutoff") {
      dataset_out |>
        dplyr::mutate(
          TICxIT = .data$tic * .data$it.ms,
          TICxIT_mean = mean(.data$TICxIT),
          is_new_outlier =
            .data$TICxIT < 1/agc_fold_cutoff * .data$TICxIT_mean |
            .data$TICxIT > agc_fold_cutoff * .data$TICxIT_mean
        ) |>
        dplyr::select(-"TICxIT", -"TICxIT_mean")
    } else if (method == "agc_sd_cutoff") {
      dataset_out |>
        dplyr::mutate(
          TICxIT = .data$tic * .data$it.ms,
          TICxIT_median = median(.data$TICxIT),
          TICxIT_sd = sd(.data$TICxIT),
          is_new_outlier =
            .data$TICxIT < (.data$TICxIT_median - agc_sd_cutoff * .data$TICxIT_sd) |
            .data$TICxIT > (.data$TICxIT_median + agc_sd_cutoff * .data$TICxIT_sd)
        ) |>
        dplyr::select(-"TICxIT", -"TICxIT_median", -"TICxIT_sd")
    },
    "something went wrong flagging outliers: "
  )
  
  # dataset outlier and type update
  if (!"is_outlier" %in% names(dataset_out)) {
    dataset_out$is_outlier <- FALSE
    dataset_out$outlier_type <- NA_character_
  }
  dataset_out <- dataset_out |>
    dplyr::mutate(
      is_outlier = ifelse(.data$is_new_outlier, TRUE, .data$is_outlier),
      outlier_type = ifelse(.data$is_new_outlier, method_type, .data$outlier_type)
    ) |>
    dplyr::select(-"is_new_outlier")
  
  # info
  n_scans_removed <- dataset_out |> dplyr::filter(.data$is_outlier) |> count_grouped_distinct("scan.no")
  if(n_scans_removed > 0){
    sprintf(
      "flagged %d/%d scans (%.1f%%) as outliers (use `orbi_plot_raw_data(y = tic * it.ms)` to visualize them) across all data groups",
      n_scans_removed, n_scans, n_scans_removed/n_scans * 100) |>
      message_finish(start_time = start_time)
  } else {
    message_finish("confirmed there are no outliers based on this method", start_time = start_time)
  }
  
  # return with restored groupings from original dataset
  return(dataset_out |> group_by_same_groups(dataset))
}

#define basepeak including nio
orbi_define_basepeak_nio <- function(dataset, basepeak_def) {
  
  # safety checks
  stopifnot(
    "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
    "`dataset` requires columns `filename`, `compound`, `scan.no`, `isotopocule`, and `ions.incremental`" =
      all(c("filename", "compound", "scan.no", "isotopocule", "ions.incremental") %in% names(dataset)),
    "`basepeak_def` needs to be a single text value identifying the isotopocule to use as the basepeak" =
      !missing(basepeak_def) && rlang::is_scalar_character(basepeak_def),
    "`basepeak_def` is not an isotopocule in the dataset" = basepeak_def %in% levels(dataset$isotopocule)
  )
  
  # ensure factors
  dataset <- dataset |> factorize_dataset("isotopocule")
  
  # info message
  start_time <-
    sprintf(
      "orbi_define_basepeak_nio() is setting the '%s' isotopocule as the ratio denominator... ",
      basepeak_def) |>
    message_start()
  
  # identify `basepeak` for each scan
  df.out <-
    try_catch_all({
      df.out <- dataset |>
        dplyr::mutate(basepeak = !!basepeak_def,
                      nio = ((.data$intensity/.data$peakNoise)*(4.4 /1))*((120000/.data$resolution)^0.5)*(.data$microscans^0.5)) |>
        dplyr::group_by(.data$filename, .data$compound, .data$scan.no)
      # add basepeak ions
      if ("is_satellite_peak" %in% names(dataset)) {
        # with satellite peak defined
        df.out <- df.out |>
          dplyr::mutate(
            basepeak_ions = .data$ions.incremental[.data$isotopocule == !!basepeak_def & !.data$is_satellite_peak],
            basepeak_nio = .data$nio[.data$isotopocule == !!basepeak_def & !.data$is_satellite_peak]
          )
      } else {
        # without satellite peak defined
        df.out <- df.out |>
          dplyr::mutate(
            basepeak_ions = .data$ions.incremental[.data$isotopocule == !!basepeak_def],
            basepeak_nio = .data$nio[.data$isotopocule == !!basepeak_def]
          )
      }
      # return
      df.out |>
        dplyr::ungroup()
    },
    # error catching
    function(p) {
      # analyze scans (is slow so only done if error)
      df_summary <-
        dataset |>
        dplyr::summarize(
          n_bp = sum(.data$isotopocule == !!basepeak_def),
          .by = c("filename", "compound", "scan.no")
        ) |>
        dplyr::summarize(
          n_scans = dplyr::n(),
          n_too_few = sum(.data$n_bp == 0L),
          n_too_many = sum(.data$n_bp > 1L),
          .by = c("filename", "compound")
        )
      
      # check abundances
      too_many_bps <- df_summary |>
        dplyr::filter(.data$n_too_many > 0)
      too_few_bps <- df_summary |>
        dplyr::filter(.data$n_too_few > 0)
      
      if (nrow(too_many_bps) > 0 && !"is_satellite_peak" %in% names(dataset)) {
        # too many base peaks
        sprintf("the %s isotopocule exists multiple times in some scans, make sure to run orbi_flag_satellite_peaks() first",
                basepeak_def) |>
          abort()
      } else if (nrow(too_few_bps) > 0) {
        info <-
          too_few_bps |>
          dplyr::mutate(
            label = sprintf("basepeak '%s' is missing in %d scans (%.1f%%) of compound '%s' in file '%s'",
                            basepeak_def, .data$n_too_few, .data$n_too_few/.data$n_scans * 100, .data$compound, .data$filename)
          ) |>
          dplyr::pull(.data$label)
        
        sprintf("the '%s' isotopocule does not exist in some scans, consider using `orbi_filter_isox()` to focus on specific file(s) and/or compound(s): \n - %s", basepeak_def,
                paste(info, collapse = "\n - ")) |>
          abort()
      } else {
        # some other error
        abort("something went wrong identifying the base peak for each scan", parent = p)
      }
    }
    )
  
  # remove basepeak from isotopocule list
  df.out <-
    try_catch_all(
      df.out |> dplyr::filter(.data$isotopocule != basepeak_def) |> droplevels(),
      "something went wrong removing the base peak isotopocule"
    )
  
  # calculate ratios
  df.out <-
    try_catch_all(
      df.out |>
        dplyr::mutate(
          ratio = .data$ions.incremental / .data$basepeak_ions,
          nio.ratio = .data$nio / .data$basepeak_nio,
          .after = "ions.incremental"
        ),
      "something went wrong calculating ratios: "
    )
  
  # info message
  sprintf(
    "set base peak and calculated %d ratios for %d isotopocules/base peak (%s)",
    nrow(df.out),
    length(levels(df.out$isotopocule)),
    paste(levels(df.out$isotopocule), collapse = ", ")
  ) |> message_finish(start_time = start_time)
  
  # return
  return(df.out)
}

# Returns a results summary table for nio values
orbi_summarize_nio_results <- function(
    dataset,
    ratio_method = c("mean", "sum", "median", "geometric_mean", "slope", "weighted_sum"),
    .by = c("block", "sample_name", "segment", "data_group", "data_type", "injection"),
    include_flagged_data = FALSE, include_unused_data = FALSE) {
  
  # basic checks
  if (missing(dataset))
    stop("no input for dataset supplied", call. = TRUE)
  
  if (is.data.frame(dataset) == FALSE)
    stop("dataset must be a data frame",  call. = TRUE)
  
  if (missing(ratio_method))
    stop("no input for ratio_method supplied", call. = TRUE)
  
  # make sure the ratio method argument is valid
  ratio_method <- rlang::arg_match(ratio_method)
  
  # check that required data columns are present
  base_group_cols <- c("filename", "compound", "basepeak", "isotopocule")
  req_cols <- c(base_group_cols, "time.min", "ions.incremental", "basepeak_ions","tic")
  if (length(missing_cols <- setdiff(req_cols, names(dataset))) > 0) {
    paste0("Missing required column(s): ",
           paste(missing_cols, collapse = ", ")) |>
      stop(call. = FALSE)
  }
  
  # filter flagged data
  n_all <- nrow(dataset)
  dataset_wo_flagged <- dataset |>
    filter(is_satellite_peak == FALSE,
           is_outlier == FALSE,
           is_weak_isotopocule == FALSE)
  n_flagged <- n_all - nrow(dataset_wo_flagged)
  if (!include_flagged_data)
    dataset <- dataset_wo_flagged
  
  # filter unused data
  n_unused <- 0L
  if ("data_type" %in% names(dataset)) {
    dataset_wo_unused <- dataset |> 
      dplyr::filter(.data$data_type == setting("data_type_data"))
    n_unused <- n_all - n_flagged - nrow(dataset_wo_unused)
    if (!include_unused_data)
      dataset <- dataset_wo_unused
  }
  
  # set basic groupings (use .add so prior group_by groupings are also preserved)
  df.group <- dataset |> dplyr::group_by(!!!lapply(base_group_cols, rlang::sym), .add = TRUE)
  
  # add additional groupings
  add_groups <- c()
  if (!missing(.by)) {
    # user defined
    by_quo <- rlang::enquo(.by)
    add_groups <- tidyselect::eval_select(expr = by_quo, data = dataset) |> names()
  } else {
    # default .by - only use the columns that actually exist in the data
    add_groups <- intersect(.by, names(dataset))
  }
  if (length(add_groups) > 0) {
    df.group <- df.group |>
      dplyr::group_by(!!!lapply(add_groups, rlang::sym), .add = TRUE)
  }
  
  # info message
  start_time <-
    sprintf(
      "orbi_summarize_nio_results() is grouping the data by %s and summarizing ratios from %d peaks (%s %d flagged peaks; %s %d unused peaks) using the '%s' method...",
      sprintf("'%s'", dplyr::group_vars(df.group)) |> paste(collapse = ", "),
      if(include_flagged_data && include_unused_data) n_all 
      else if (include_flagged_data) n_all - n_unused
      else if (include_unused_data) n_all - n_flagged
      else n_all - n_unused - n_flagged, 
      if(include_flagged_data) "including" else "excluded", n_flagged,
      if(include_unused_data) "including" else "excluded", n_unused,
      ratio_method
    ) |>
    message_start()
  
  # calculations
  df.stat <-
    try_catch_all({
      # run calculations
      df.group |>
        
        dplyr::summarize(
          
          # scan information
          start_scan.no = min(.data$scan.no),
          end_scan.no = max(.data$scan.no),
          
          # time information
          start_time.min = min(.data$time.min),
          mean_time.min = mean(.data$time.min),
          end_time.min = max(.data$time.min),
          
          #scan info
          number_of_scans = length(.data$ions.incremental / .data$basepeak_ions),
          coverage = number_of_scans/max(.data$scan.no)*100, 
          
          tic_m= mean(tic),
          tic_sd= sd(tic),
          tic_sem=sd(tic)/sqrt(length(tic)),
          
          #isotopocule info
          ions.inc.m=mean(ions.incremental),
          ions.inc.sd=sd(ions.incremental),
          ions.inc.sem=sd(ions.incremental)/sqrt(length(ions.incremental)),
          basepeak.ions.m=mean(basepeak_ions),
          basepeak.ions.sd=sd(basepeak_ions),
          basepeak.ions.sem=sd(basepeak_ions)/sqrt(length(basepeak_ions)),
          nio.m=mean(nio),
          nio.sd=sd(nio),
          nio.sem=sd(nio)/sqrt(length(nio)),
          basepeak.nio.m=mean(basepeak_nio),
          basepeak.nio.sd=sd(basepeak_nio),
          basepeak.nio.sem=sd(basepeak_nio)/sqrt(length(basepeak_nio)),
          
          # ratio calculation
          ratio = orbi_calculate_summarized_ratio(
            .data$ions.incremental,
            .data$basepeak_ions,
            ratio_method = !!ratio_method
          ),
          nio.ratio= orbi_calculate_summarized_ratio(
            .data$nio,
            .data$basepeak_nio,
            ratio_method = !!ratio_method
          ),
          # standard error calculation 
          ratio_sem = calculate_ratios_sem(
            ratios = .data$ions.incremental / .data$basepeak_ions),
          nio.ratio_sem = calculate_ratios_sem(
            ratios = .data$nio / .data$basepeak_nio),
          
          # data quality measures
          minutes_to_1e6_ions = (1E6 / sum(.data$ions.incremental)) *
            (max(.data$time.min) - min(.data$time.min)),
          minutes_to_1e6_ions_nio = (1E6 / sum(.data$nio)) *
            (max(.data$time.min) - min(.data$time.min)),
          
          shot_noise_permil =
            1000 * (sqrt(
              (sum(.data$ions.incremental) + sum(.data$basepeak_ions)) /
                (sum(.data$ions.incremental) * sum(.data$basepeak_ions)))),
          shot_noise_permil_nio =
            1000 * (sqrt(
              (sum(.data$nio) + sum(.data$basepeak_nio)) /
                (sum(.data$nio) * sum(.data$basepeak_nio)))),
          
          .groups = "drop") |>
        
        dplyr::mutate(ratio_relative_sem_permil = 1000 * (.data$ratio_sem / .data$ratio),
                      ratio_relative_sem_permil_nio = 1000 * (.data$nio.ratio_sem / .data$nio.ratio))   |>
        
        # round values for output
        dplyr::mutate(
          ions.inc.m=round(.data$ions.inc.m, 8),
          ions.inc.sd=round(.data$ions.inc.sd, 8),
          ions.inc.sem=round(.data$ions.inc.sem, 8),
          basepeak.ions.m=round(.data$basepeak.ions.m, 8),
          basepeak.ions.sd=round(.data$basepeak.ions.sd, 8),
          basepeak.ions.sem=round(.data$basepeak.ions.sem, 8),
          nio.m=round(.data$nio.m, 8),
          nio.sd=round(.data$nio.sd, 8),
          nio.sem=round(.data$nio.sem, 8),
          basepeak.nio.m=round(.data$basepeak.nio.m, 8),
          basepeak.nio.sd=round(.data$basepeak.nio.sd, 8),
          basepeak.nio.sem=round(.data$basepeak.nio.sem, 8),
          
          ratio = round(.data$ratio, 8),
          ratio_sem = round(.data$ratio_sem, 8),
          ratio_relative_sem_permil = round(.data$ratio_relative_sem_permil, 3),
          shot_noise_permil = round(.data$shot_noise_permil, 3),
          SNL_limit= ifelse(.data$ratio_relative_sem_permil<=2*shot_noise_permil, "pass","fail"),
          minutes_to_1e6_ions = round(.data$minutes_to_1e6_ions, 2),
          
          nio.ratio = round(.data$nio.ratio, 8),
          nio.ratio_sem = round(.data$nio.ratio_sem, 8),
          ratio_relative_sem_permil_nio = round(.data$ratio_relative_sem_permil_nio, 3),
          shot_noise_permil_nio = round(.data$shot_noise_permil_nio, 3),
          SNL_limit_nio= ifelse(.data$ratio_relative_sem_permil_nio<=2*shot_noise_permil_nio, "pass","fail"),
          minutes_to_1e6_ions_nio = round(.data$minutes_to_1e6_ions_nio, 2),
          
          tic_m= round(.data$tic_m, 10),
          tic_sd= round(.data$tic_sd, 8),
          tic_sem=round(.data$tic_sem, 8),
          
          coverage=round(.data$coverage, 2),
        )  |>
        # sort table by the grouping variables
        dplyr::arrange(!!!lapply(dplyr::group_vars(df.group), rlang::sym)) 
    },
    "something went wrong summarizing the results: "
    )
  
  # info
  message_finish("completed", start_time = start_time)
  
  # return
  return(df.stat)
  
}

#define nio and add isotopocule counts in seperate column
#run after basepeak filter for any scan that is missing unsubsituted ion 
orbi_ro_data <- function(dataset, basepeak_unsub) {
  
  # safety checks
  stopifnot(
    "need a `dataset` data frame" = !missing(dataset) && is.data.frame(dataset),
    "`dataset` requires columns `filename`, `compound`, `scan.no`, `isotopocule`, and `ions.incremental`" =
      all(c("filename", "compound", "scan.no", "isotopocule", "ions.incremental") %in% names(dataset)),
    "`basepeak_unsub` needs to be a single text value identifying the isotopocule to use as the basepeak" =
      !missing(basepeak_unsub) && rlang::is_scalar_character(basepeak_unsub),
    "`basepeak_unsub` is not an isotopocule in the dataset" = basepeak_unsub %in% levels(dataset$isotopocule)
  )
  
  # ensure factors
  dataset <- dataset |> factorize_dataset("isotopocule")

  df.out <- dataset |>
    dplyr::mutate(basepeak = !!basepeak_unsub,
                  nio = ((intensity/peakNoise)*(4.4 /1))*((120000/resolution)^0.5)*(microscans^0.5)) |>
    dplyr::group_by(filename, compound, scan.no)|>
    dplyr::mutate(
      U0_ions = ions.incremental[isotopocule == !!basepeak_unsub],
      U0_nio = nio[isotopocule == !!basepeak_unsub],
      H2_ions = ifelse(any(isotopocule == "2H"),ions.incremental[isotopocule == "2H"],NA), #b/c 0 row values
      H2_nio = ifelse(any(isotopocule == "2H"),nio[isotopocule == "2H" ],NA),
      O17_ions = ifelse(any(isotopocule == "17O"),ions.incremental[isotopocule == "17O"],NA),
      O17_nio = ifelse(any(isotopocule == "17O"),nio[isotopocule == "17O"],NA)
    )|>
    #remove isotopocules that are already defined in columns
    dplyr::ungroup()|> 
    dplyr::filter(isotopocule != basepeak_unsub) |>
    dplyr::filter(isotopocule != "2H")|> 
    dplyr::filter(isotopocule != "17O")|> 
    dplyr::mutate(
      ratio.13C.0U = ions.incremental / U0_ions,
      ratio.nio.13C.0U = nio / U0_nio)|>
    dplyr::rename(C13_ions=ions.incremental,
                  C13_nio=nio)
  # return
  return(df.out)
}

# Returns a results summary table for ro value table
#only works after orbi_ro_data () and renaming ions.incrmental column and nio column to C13_ions, C13_nio 
orbi_ro_summarize_results <- function(dataset,ratio_method) {
  
  base_group_cols <- c("filename", "compound", "basepeak")
  
  # set basic groupings (use .add so prior group_by groupings are also preserved)
  df.group <- dataset |> dplyr::group_by(!!!lapply(base_group_cols, rlang::sym), .add = TRUE)
  
  # calculations
  df.stat <-
    try_catch_all({
      # run calculations
      df.group |>
        
        dplyr::summarize(
          #tic data
          tic_m= mean(tic),
          tic_sd= sd(tic),
          tic_sem=sd(tic)/sqrt(length(tic)),
          
          #isotopocule data
          m.U0_ions=mean(U0_ions),
          m.C13_ions=mean(C13_ions),
          m.H2_ions=mean(H2_ions),
          m.O17_ions=mean(O17_ions),
          
          sd.U0_ions=sd(U0_ions),
          sd.C13_ions=sd(C13_ions),
          sd.H2_ions=sd(H2_ions),
          sd.O17_ions=sd(O17_ions),
          
          sem.U0_ions=sd(U0_ions)/sqrt(length(U0_ions)),
          sem.C13_ions=sd(C13_ions)/sqrt(length(C13_ions)),
          sem.H2_ions=sd(H2_ions)/sqrt(length(H2_ions)),
          sem.O17_ions=sd(O17_ions)/sqrt(length(O17_ions)),
          
          rse.U0_ions=1000 * (sem.U0_ions/m.U0_ions),
          rse.C13_ions=1000 * (sem.C13_ions/m.C13_ions),
          rse.H2_ions=1000 * (sem.H2_ions/m.H2_ions),
          rse.O17_ions=1000 * (sem.O17_ions/m.O17_ions),
          
          m.U0_nio=mean(U0_nio),
          m.C13_nio=mean(C13_nio),
          m.H2_nio=mean(H2_nio),
          m.O17_nio=mean(O17_nio),
          
          sd.U0_nio=sd(U0_nio),
          sd.C13_nio=sd(C13_nio),
          sd.H2_nio=sd(H2_nio),
          sd.O17_nio=sd(O17_nio),
          
          sem.U0_nio=sd(U0_nio)/sqrt(length(U0_nio)),
          sem.C13_nio=sd(C13_nio)/sqrt(length(C13_nio)),
          sem.H2_nio=sd(H2_nio)/sqrt(length(H2_nio)),
          sem.O17_nio=sd(O17_nio)/sqrt(length(O17_nio)),
          
          rse.U0_nio=1000 * (sem.U0_nio/m.U0_nio),
          rse.C13_nio=1000 * (sem.C13_nio/m.C13_nio),
          rse.H2_nio=1000 * (sem.H2_nio/m.H2_nio),
          rse.O17_nio=1000 * (sem.O17_nio/m.O17_nio),
          
          # ratio calculation
          ratio = orbi_calculate_summarized_ratio(
            .data$C13_ions,
            .data$U0_ions,
            ratio_method = !!ratio_method),
          nio.ratio= orbi_calculate_summarized_ratio(
            .data$C13_nio,
            .data$U0_nio,
            ratio_method = !!ratio_method),
          # standard error calculation 
          ratio_sem = calculate_ratios_sem(
            ratios = .data$C13_ions / .data$U0_ions),
          nio.ratio_sem = calculate_ratios_sem(
            ratios = .data$C13_nio / .data$U0_nio),
          
          #ro calculations
          ro13C.0U.ions=m.C13_ions/
            sum(m.C13_ions,m.U0_ions),
          prop.err.ro13C.0U.ions=sqrt(abs(((m.U0_ions/(m.C13_ions+m.U0_ions)^2)*sem.C13_ions)^2+((-(m.C13_ions/(m.C13_ions+m.U0_ions)^2))))),
          ro13C.0U.nio=m.C13_nio/
            sum(m.C13_nio,m.U0_nio),
          prop.err.ro13C.0U.nio=sqrt(abs(((m.U0_nio/(m.C13_nio+m.U0_nio)^2)*sem.C13_nio)^2+((-(m.C13_nio/(m.C13_nio+m.U0_nio)^2))))),
          
          .groups = "drop") |>
        # sort table by the grouping variables
        dplyr::arrange(!!!lapply(dplyr::group_vars(df.group), rlang::sym))
      
    },
    "something went wrong summarizing the results: "
    )
  
  # return
  return(df.stat)
  
}

orbi_ro_output <- function(dataset) {
  
  base_group_cols <- c("filename", "compound", "basepeak")
  
  # set basic groupings (use .add so prior group_by groupings are also preserved)
  df.group <- dataset |> dplyr::group_by(!!!lapply(base_group_cols, rlang::sym), .add = TRUE)
  
  # calculations
  df.stat <-
    try_catch_all({
      # run calculations
      df.group |>
        
        dplyr::summarize(
          ro13C.0U.ions=m.C13_ions/
            sum(m.C13_ions,m.U0_ions),
          prop.err.ro13C.0U.ions=sqrt(abs(((m.U0_ions/(m.C13_ions+m.U0_ions)^2)*sem.C13_ions)^2+((-(m.C13_ions/(m.C13_ions+m.U0_ions)^2))))),
          ro13C.0U.nio=m.C13_nio/
            sum(m.C13_nio,m.U0_nio),
          prop.err.ro13C.0U.nio=sqrt(abs(((m.U0_nio/(m.C13_nio+m.U0_nio)^2)*sem.C13_nio)^2+((-(m.C13_nio/(m.C13_nio+m.U0_nio)^2))))),
          .groups = "drop") |>
        # sort table by the grouping variables
        dplyr::arrange(!!!lapply(dplyr::group_vars(df.group), rlang::sym))
      
    },
    "something went wrong summarizing the results: "
    )
  
  # return
  return(df.stat)
  
}

### fractional abundance
f.abundance <- function(x){
  data= x %>%
    mutate(element = str_extract(isotopocule, "(?<=\\d)\\p{L}+"))%>%
    group_by(filename,compound, isotopocule) %>%
    mutate(n.total.scans =n())%>%
    ungroup()
  data_out = data[order(data$filename,data$compound,data$isotopocule,data$scan.no),]#sort data
  #calculate fractional abundances 
  data_out=data_out %>%
    #calculate how much percent of scans have 2H relative to 0U per compound
    group_by(filename, compound) %>% #count occurrences of ilog per compound/fragment group
    mutate(rel.ab.2H=length(which(isotopocule == "2H"))/length(unique(scan.no))*100)%>%
    ungroup()%>%
    #calculate Nio
    mutate(Nio = ((intensity/peakNoise)*(4.4 /1))*((120000/resolution)^0.5)*(microscans^0.5))%>% #append column with calculated Nio values (Eiler 2017) CN=4.4 ref.resolution=120000, z=1)
    #create simulated dataset
    #Identifies missing isotopocules of interest per scan. Needed to add Nio value wherever 13C or 0U missing
    complete(nesting(filename,compound,scan.no),
             isotopocule,
             fill = list(Nio=0),
             explicit = FALSE)%>%
    mutate(simulated = is.na(intensity),
           element = case_when(simulated == "TRUE" & isotopocule == "0U" ~ "U",
                               simulated == "TRUE" & isotopocule == "13C" ~ "C",
                               simulated == "TRUE" & isotopocule == "2H" ~ "H",
                               simulated == "TRUE" & isotopocule == "17O" ~ "O",
                               TRUE ~ element))%>%
    filter(!grepl("O",compound) & element != "O"| grepl("O",compound))%>%#removing O rows for fragments that do not have O isotopocule
    mutate(Nio = ifelse(simulated == "TRUE" & element == "U",1,Nio))%>%
    mutate(n.H=case_when(compound=="C5H7"~ 7,compound=="C5H9"|compound=="C5H9O"|compound=="C6H9"~ 9,
                         compound=="C5H11"|compound=="C6H11"|compound=="C6H11O"|compound=="C7H11"~ 11,
                         compound=="C6H13"|compound=="C7H13"|compound=="C7H13O"|compound=="C8H13"~ 13,
                         compound=="C7H15"|compound=="C8H15"|compound=="C8H15O"|compound=="C9H15"~ 15,
                         compound=="C8H17"|compound=="C9H17"|compound=="C9H17O"|compound=="C10H17"~ 17,
                         compound=="C9H19"|compound=="C10H19"|compound=="C11H19"~ 19,
                         compound=="C11H21"|compound=="C12H21"~ 21,
                         compound=="C12H23"~ 23),
           n.C=case_when(compound=="C5H7"|compound=="C5H9"|compound=="C5H9O"|compound=="C5H11"~ 5,
                         compound=="C6H9"|compound=="C6H11"|compound=="C6H13"|compound=="C6H11O"~ 6,
                         compound=="C7H11"|compound=="C7H13"|compound=="C7H15"|compound=="C7H13O"~ 7,
                         compound=="C8H17"|compound=="C8H13"|compound=="C8H15"|compound=="C8H15O"~ 8,
                         compound=="C9H15"|compound=="C9H17"|compound=="C9H19"|compound=="C9H17O"~ 9,
                         compound=="C10H17"|compound=="C10H19"~ 10,
                         compound=="C11H19"|compound=="C11H21"~ 11,
                         compound=="C12H21"|compound=="C12H23"~ 12))%>%
    #calculate fa denom not including 2H
    group_by(filename,compound,scan.no) %>%
    mutate(fa_denom_sum_Nio_noH = ifelse(isotopocule != "2H", sum(Nio), NA),
           fa_Nio_noH= Nio/fa_denom_sum_Nio_noH, NA)%>%
    ungroup()%>%
    #calculate fa including H
    group_by(filename,compound,scan.no)%>%
    mutate(fa_denom_sum_Nio= sum(Nio),# calculate fractional abundance denominator based on sum of Nio value of all ilog of interest in scan
           fa_Nio=  Nio/fa_denom_sum_Nio)%>%
    ungroup()%>%
    #calculate delta value
    group_by(filename,compound,scan.no)%>%
    mutate(delta_C_R=abs((deltaC/1000+1)*0.01118)/(1+(deltaC/1000+1)*0.01118),
           delta_H_R=abs((deltaH/1000+1)*0.000156)/(1+(deltaH/1000+1)*0.000156),
           delta_H_R_OH=abs((deltaOH/1000+1)*0.000156)/(1+(deltaOH/1000+1)*0.000156),
           delta_H_R_APCI=abs((deltaAPCI/1000+1)*0.000156)/(1+(deltaAPCI/1000+1)*0.000156),
           delta_H_R_O=abs((deltaO/1000+1)*0.000381)/(1+(deltaO/1000+1)*0.000381),
           delta_H_R_blank=abs((0/1000+1)*0.002005)/(1+(0/1000+1)*0.002005))%>%
    mutate(corr_factor= (((1-delta_C_R)^20)*
                           ((1-delta_H_R)^37)*
                           (1-delta_H_R_OH)*
                           (1-delta_H_R_APCI)*
                           (1-delta_H_R_O-delta_H_R_blank)))%>%
    mutate(corr_factor2= 20*abs(corr_factor/(1-delta_C_R)*delta_C_R)+
             37*abs(corr_factor/(1-delta_H_R)*delta_H_R)+
             abs(corr_factor/(1-delta_H_R_OH)*delta_H_R_OH)+
             abs(corr_factor/(1-delta_H_R_APCI)*delta_H_R_APCI)+
             abs(corr_factor/(1-delta_H_R_O-delta_H_R_blank)*delta_H_R_O))%>%
    mutate(U_M1=corr_factor2/corr_factor,
           H2_corr=(abs(corr_factor/(1-delta_H_R)*delta_H_R)/corr_factor2))%>%
    select(filename,ID, date, compound, isotopocule, element, n.H, n.C, scan.no, time.min,m_z, n.total.scans, tic, it.ms,intensity,
           ions.incremental, Nio, fa_Nio, fa_Nio_noH, U_M1, H2_corr,rel.ab.2H, simulated)%>%
    arrange(filename,isotopocule,n.C,n.H)
  return(data_out)
}
