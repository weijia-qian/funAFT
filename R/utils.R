# ---- helpers ----
`%||%` <- function(a, b) if (is.null(a)) b else a

.is_df <- function(obj) inherits(obj, "data.frame")

# Fetch a single vector by name from data or env
.get_vec <- function(nm, data, env, required_n = NULL, label = nm) {
  val <- NULL
  if (!is.null(data) && nm %in% names(data)) val <- data[[nm]]
  if (is.null(val) && exists(nm, envir = env, inherits = TRUE)) val <- get(nm, envir = env)
  if (is.null(val)) stop("Could not find variable `", nm, "` in `data` or the environment.")

  # If it's a 1D vector, accept; if it's a 2D object with 1 column, drop to vector
  if (is.matrix(val) || .is_df(val)) {
    if (ncol(as.matrix(val)) != 1L) {
      stop("`", label, "` must be a vector; got a matrix/data.frame with ", ncol(as.matrix(val)), " columns.")
    }
    val <- as.vector(val[, 1])
  } else {
    val <- as.vector(val)
  }
  if (!is.null(required_n) && length(val) != required_n) {
    stop("Length mismatch for `", label, "`: expected ", required_n, ", got ", length(val), ".")
  }
  val
}

# Bind one or many columns named in `nms`, where each item can be vector or matrix/data.frame in data or env
.bind_cols <- function(nms, data, env, required_n = NULL, label = "X") {
  pieces <- list()
  for (nm in nms) {
    val <- NULL
    if (!is.null(data) && nm %in% names(data)) val <- data[[nm]]
    if (is.null(val) && exists(nm, envir = env, inherits = TRUE)) val <- get(nm, envir = env)
    if (is.null(val)) stop("Could not find `", nm, "` for `", label, "` in `data` or the environment.")

    mat <- NULL
    if (is.vector(val)) {
      mat <- matrix(as.vector(val), ncol = 1L)
      colnames(mat) <- nm
    } else if (is.matrix(val)) {
      mat <- val
      if (is.null(colnames(mat))) colnames(mat) <- paste0(nm, "_", seq_len(ncol(mat)))
    } else if (.is_df(val)) {
      mat <- as.matrix(val)
      if (is.null(colnames(mat))) colnames(mat) <- paste0(nm, "_", seq_len(ncol(mat)))
    } else {
      stop("Object `", nm, "` for `", label, "` must be a vector, matrix, or data.frame.")
    }

    if (!is.null(required_n) && nrow(mat) != required_n) {
      stop("Row mismatch for `", nm, "` in `", label, "`: expected ", required_n, ", got ", nrow(mat), ".")
    }
    pieces[[nm]] <- mat
  }
  out <- do.call(cbind, pieces)
  # If duplicate colnames emerged, make.unique
  colnames(out) <- make.unique(colnames(out))
  out
}

# Regex discovery for X colnames across data or env
.discover_x_by_regex <- function(pattern, data, env) {
  col_hits <- character(0)
  if (!is.null(data)) {
    col_hits <- grep(pattern, names(data), value = TRUE)
  }
  # If none in data, try environment object names
  if (!length(col_hits)) {
    cand <- ls(envir = env, all.names = TRUE)
    obj_hits <- grep(pattern, cand, value = TRUE)
    # Keep only objects that are vector/matrix/df (exclude functions, etc.)
    keep <- vapply(obj_hits, function(nm) {
      obj <- get(nm, envir = env)
      is.vector(obj) || is.matrix(obj) || .is_df(obj)
    }, logical(1))
    obj_hits[keep]
  } else {
    col_hits
  }
}
