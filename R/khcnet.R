#' KH Coder-style Co-occurrence Network (quanteda + tidygraph + ggraph)
#'
#' Build and plot a KH Coder-like co-occurrence network from a quanteda dfm.
#' - Community detection (fast_greedy / walktrap / edge_betweenness)
#' - Within-community vs between-community edges (solid vs dotted)
#' - Grayscale-print mode (node fill in greys)
#' - Optional node-size (frequency) legend
#' - Optional node fill transparency in color mode
#' - CJK font family and label color options
#'
#' @param dfm A quanteda dfm (doc x term).
#' @param community Community detection method:
#'   "fast_greedy", "walktrap", "edge_betweenness", "none".
#' @param layout Graph layout: "fr" or "kk".
#' @param seed Random seed for reproducible layout.
#' @param edges Keep top-N edges by weight. If not NULL, overrides th.
#' @param th Numeric threshold for edge weights. Used only when edges is NULL.
#' @param remove_isolates Remove isolated nodes (degree 0) before plotting.
#' @param use_freq_as_size If TRUE, node size is proportional to term frequency.
#' @param bubble_size Multiplier for node size (affects max size).
#' @param show_size_legend If TRUE, show node-size (frequency) legend.
#' @param bw_groups If TRUE, grayscale printing mode.
#' @param bw_intercommunity_linetype In bw_groups mode, inter-community edge linetype:
#'   "dotted" (default) or "solid".
#' @param fill_alpha Node fill alpha in color mode (0–1). Ignored if bw_groups = TRUE.
#' @param font_family Font family for node labels.
#' @param label_colour Color for node labels.
#' @param label_repel If TRUE, repel node labels.
#' @param label_size Numeric label size.
#'
#' @return A list with graph, tbl_graph, plot, info.
#' @export
khcnet <- function(
    dfm,
    community = c("fast_greedy", "walktrap", "edge_betweenness", "none"),
    layout = c("fr", "kk"),
    seed = 100,
    edges = 80,
    th = NULL,
    remove_isolates = TRUE,
    use_freq_as_size = TRUE,
    bubble_size = 100,
    show_size_legend = FALSE,
    bw_groups = FALSE,
    bw_intercommunity_linetype = c("dotted", "solid"),
    fill_alpha = 0.8,
    font_family = "sans",
    label_colour = "black",
    label_repel = TRUE,
    label_size = 3
){

  community <- match.arg(community)
  layout <- match.arg(layout)
  bw_intercommunity_linetype <- match.arg(bw_intercommunity_linetype)

  if (!requireNamespace("quanteda", quietly = TRUE)) stop("Package 'quanteda' is required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  if (!requireNamespace("tidygraph", quietly = TRUE)) stop("Package 'tidygraph' is required.")
  if (!requireNamespace("ggraph", quietly = TRUE)) stop("Package 'ggraph' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  if (!quanteda::is.dfm(dfm)) {
    stop("`dfm` must be a quanteda dfm (doc x term).")
  }

  if (!is.numeric(fill_alpha) || length(fill_alpha) != 1 || is.na(fill_alpha) ||
      fill_alpha < 0 || fill_alpha > 1) {
    stop("`fill_alpha` must be a single number in [0, 1].")
  }

  # ---- term-term co-occurrence ----
  fcm <- quanteda::fcm(dfm)
  W <- Matrix::Matrix(fcm, sparse = TRUE)
  class(W) <- "dgCMatrix"

  g0 <- igraph::graph_from_adjacency_matrix(
    W, mode = "undirected", weighted = TRUE, diag = FALSE
  )

  if (remove_isolates) {
    g0 <- igraph::induced_subgraph(g0, igraph::degree(g0) > 0)
  }

  if (igraph::ecount(g0) == 0) {
    stop("No edges to draw. Try a less sparse dfm or adjust preprocessing.")
  }

  el <- igraph::as_data_frame(g0, what = "edges")
  if (nrow(el) == 0) stop("No edges to draw (edge table is empty).")

  # ---- threshold from edges OR th ----
  if (!is.null(edges)) {
    edges <- as.integer(edges)
    if (is.na(edges) || edges < 1) stop("`edges` must be a positive integer.")

    # ★重要：edges を存在するエッジ数で上限カット
    edges <- min(edges, nrow(el))

    # probs を必ず[0,1]に収める
    probs <- 1 - edges / nrow(el)
    probs <- max(0, min(1, probs))

    th_use <- as.numeric(stats::quantile(el$weight, probs = probs, names = FALSE))
  } else {
    if (is.null(th)) stop("Either `edges` or `th` must be provided.")
    th_use <- as.numeric(th)
    if (!is.finite(th_use)) stop("`th` must be a finite numeric value.")
  }

  el2 <- el[el$weight >= th_use, , drop = FALSE]
  if (nrow(el2) == 0) stop("No edges left after filtering. Try larger `edges` or smaller `th`.")

  g <- igraph::graph_from_data_frame(
    el2,
    directed = FALSE,
    vertices = igraph::as_data_frame(g0, what = "vertices")
  )

  if (remove_isolates) {
    g <- igraph::induced_subgraph(g, igraph::degree(g) > 0)
  }

  if (igraph::ecount(g) == 0 || igraph::vcount(g) < 2) {
    stop("No drawable graph after filtering/removing isolates.")
  }

  # ---- node frequency ----
  freq_all <- quanteda::featfreq(dfm)
  freq <- as.numeric(freq_all[igraph::V(g)$name])
  freq[is.na(freq)] <- 0

  igraph::V(g)$freq <- freq
  igraph::V(g)$freq_plot <- if (use_freq_as_size) freq else 1
  igraph::V(g)$comm <- 1L

  # ---- community detection ----
  if (community == "fast_greedy") {
    com <- igraph::cluster_fast_greedy(g, weights = igraph::E(g)$weight)
    igraph::V(g)$comm <- igraph::membership(com)
  } else if (community == "walktrap") {
    com <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight)
    igraph::V(g)$comm <- igraph::membership(com)
  } else if (community == "edge_betweenness") {
    com <- igraph::cluster_edge_betweenness(g, weights = igraph::E(g)$weight)
    igraph::V(g)$comm <- igraph::membership(com)
  } else {
    igraph::V(g)$comm <- 1L
  }
  igraph::V(g)$comm <- factor(igraph::V(g)$comm)

  # ---- edge linetype (within vs between) ----
  ends <- igraph::ends(g, igraph::E(g), names = TRUE)
  idx <- setNames(seq_along(igraph::V(g)$name), igraph::V(g)$name)
  a <- idx[ends[,1]]
  b <- idx[ends[,2]]
  same_comm <- igraph::V(g)$comm[a] == igraph::V(g)$comm[b]

  if (bw_groups) {
    inter_lty <- if (bw_intercommunity_linetype == "solid") "solid" else "dotted"
    igraph::E(g)$lty <- ifelse(same_comm, "solid", inter_lty)
  } else {
    igraph::E(g)$lty <- ifelse(same_comm, "solid", "dotted")
  }

  tg <- tidygraph::as_tbl_graph(g)
  set.seed(seed)

  size_mult <- bubble_size / 100
  size_aes <- if (use_freq_as_size) "freq_plot" else 1

  p <- ggraph::ggraph(tg, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(linetype = lty),
      colour = "grey50"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = .data[[size_aes]], fill = comm),
      shape = 21,
      colour = "grey30",
      stroke = 0.8,
      alpha = if (bw_groups) 1 else fill_alpha
    ) +
    ggplot2::scale_linetype_identity() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_family = font_family)

  # ---- size legend ----
  if (show_size_legend && use_freq_as_size) {
    p <- p + ggplot2::scale_size_area(
      name = "Frequency",
      max_size = 20 * size_mult
    )
  } else {
    p <- p + ggplot2::scale_size_area(
      max_size = 20 * size_mult,
      guide = "none"
    )
  }

  # ---- labels ----
  if (label_repel && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggraph::geom_node_text(
      ggplot2::aes(label = name),
      family = font_family,
      size = label_size,
      colour = label_colour,
      repel = TRUE
    )
  } else {
    p <- p + ggraph::geom_node_text(
      ggplot2::aes(label = name),
      family = font_family,
      size = label_size,
      colour = label_colour
    )
  }

  # ---- fill scale ----
  if (bw_groups) {
    p <- p + ggplot2::scale_fill_grey(start = 0.88, end = 0.35)
  } else {
    p <- p + ggplot2::scale_fill_hue()
  }

  info <- list(
    nodes = igraph::vcount(g),
    edges = igraph::ecount(g),
    density = igraph::edge_density(g),
    threshold = th_use
  )

  list(graph = g, tbl_graph = tg, plot = p, info = info)
}
