extract_barcode_info <- function(filename) {
  # Input validation
  if (!is.character(filename) || length(filename) != 1) {
    stop("Input 'filename' must be a single character string.")
  }

  # --- Barcode Extraction ---

  # Define patterns for barcode (using PERL-compatible regex for consistency)
  barcode_patterns <- c(
    '^([ACGT]{8})_',  # Barcode at the beginning followed by underscore
    '_([ACGT]{8})_',  # Barcode flanked by two underscores
    '_([ACGT]{8})$'   # Barcode preceded by an underscore at the end
  )

  barcode <- NULL # Initialize barcode to NULL

  # Loop through barcode patterns to find the first match
  for (pattern in barcode_patterns) {
    # regexec finds the positions of the match and captured groups
    match_info <- regexec(pattern, filename, perl = TRUE)

    # Check if a match was found (match position is not -1)
    if (match_info[[1]][1] != -1) {
      # regmatches extracts the actual matched strings
      matches <- regmatches(filename, match_info)
      # The captured group is the second element
      # matches[[1]] contains c(full_match, group1, group2, ...)
      barcode <- matches[[1]][2]
      break # Stop after finding the first match
    }
  }

  # --- Barcode ID Extraction ---

  # Define patterns for barcode_id (matching the provided Python code)
  # Python code specifically looks for 'bcode_(\d+)'
  # In R, '\d' needs to be escaped: '\\d'
  barcode_id_patterns <- c(
    'bcode_(\\d+)'
    # If other patterns were intended per the Python docstring, add them here:
    # '_(\\d+)_',        # ID flanked by underscores
    # '_(\\d+)\\.[^.]+$' # ID after underscore before extension
  )

  barcode_id <- NULL # Initialize barcode_id to NULL

  # Loop through barcode_id patterns to find the first match
  for (pattern in barcode_id_patterns) {
    match_info <- regexec(pattern, filename, perl = TRUE)

    # Check if a match was found
    if (match_info[[1]][1] != -1) {
      matches <- regmatches(filename, match_info)
      # Extract the captured group (the digits string)
      barcode_id_str <- matches[[1]][2]

      # Convert the extracted string to an integer
      # Use suppressWarnings to handle potential non-numeric captures gracefully (though \d+ ensures digits)
      converted_id <- suppressWarnings(as.integer(barcode_id_str))

      # Assign if conversion was successful (not NA)
      if (!is.na(converted_id)) {
          barcode_id <- converted_id
      }
      # If conversion results in NA, barcode_id remains NULL

      break # Stop after finding the first match
    }
  }
    # Return results as a named list (R's equivalent to a tuple/dictionary)
  return(list(barcode = barcode, barcode_id = barcode_id))
}
