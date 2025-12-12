#' KH Coder-style Co-occurrence Network (quanteda + tidygraph + ggraph)
#'
#' Build and plot a KH Coder-like co-occurrence network from a quanteda dfm.
#' - Community detection (fast_greedy / walktrap / edge_betweenness)
#' - Within-community vs between-community edges (solid vs dotted) [KH Coder rule]
#' - Color mode vs grayscale printing mode (node fill)
#' - Option to control node fill transparency in color mode
#' - CJK font family and label color options
#' - Option to show size legend for node frequency
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
#' @param bw_groups If TRUE, grayscale printing mode for node fill.
#' @param fill_alpha In color mode (bw_groups = FALSE), node fill alpha in [0, 1].
#'   Ignored in grayscale mode.
#' @param show_size_legend If TRUE, show size legend for node frequency when use_freq_as_size = TRUE.
#' @param font_family Font family for node labels (e.g., "Noto Sans CJK JP", "Noto Sans CJK KR").
#' @param label_colour Color for node labels (e.g., "black", "gray20", "#333333").
#' @param label_repel If TRUE, repel node labels (requires ggrepel; falls back if missing).
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
    bw_groups = FALSE,
    fill_alpha = 0.75,
    show_size_legend = FALSE,
    font_family = "sans",
    label_colour = "black",
    label_repel = TRUE,
    label_size = 3
){
  community <- match.arg(community)
  layout <- match.arg(layout)

  if (!requireNamespace("quanteda", quietly = TRUE)) stop("Package 'quanteda' is required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  if (!requireNamespace("tidygraph", quietly = TRUE)) stop("Package 'tidygraph' is required.")
  if (!requireNamespace("ggraph", quietly = TRUE)) stop("Package 'ggraph' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  if (!quanteda::is.dfm(dfm)) {
    stop("`dfm` must be a quanteda dfm object (doc x term).")
  }
  if (is.null(dim(dfm)) || length(dim(dfm)) != 2) {
    stop("`dfm` must be 2D (doc x term). If you subset, use drop = FALSE.")
  }

  # ---- (A) term-term co-occurrence matrix (fcm) ----
  fcm <- quanteda::fcm(dfm)  # feature x feature

  # Convert fcm to a plain sparse matrix (drop quanteda classes safely)
  W <- Matrix::Matrix(fcm, sparse = TRUE)
  class(W) <- "dgCMatrix"

  # ---- (B) igraph from adjacency ----
  g0 <- igraph::graph_from_adjacency_matrix(
    W,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  if (remove_isolates) {
    keep <- igraph::degree(g0) > 0
    g0 <- igraph::induced_subgraph(g0, vids = which(keep))
  }

  if (igraph::vcount(g0) < 2 || igraph::ecount(g0) == 0) {
    stop("No edges to draw. Try lowering thresholds or check dfm content.")
  }

  # ---- (C) edge filtering: top-N edges OR threshold th ----
  el <- igraph::as_data_frame(g0, what = "edges")
  if (nrow(el) == 0) stop("No edges to draw (graph has zero edges).")

  if (!is.null(edges)) {
    edges <- as.integer(edges)
    if (is.na(edges) || edges < 1) stop("`edges` must be a positive integer.")
    edges <- min(edges, nrow(el))
    # IMPORTANT: guard against probs outside [0,1]
    probs <- 1 - edges / nrow(el)
    probs <- max(min(probs, 1), 0)
    th_use <- as.numeric(stats::quantile(el$weight, probs = probs, names = FALSE))
  } else {
    if (is.null(th)) stop("Either `edges` or `th` must be provided.")
    th_use <- as.numeric(th)
  }

  el2 <- el[el$weight >= th_use, , drop = FALSE]
  if (nrow(el2) == 0) stop("No edges left after filtering. Try larger `edges` or smaller `th`.")

  g <- igraph::graph_from_data_frame(
    el2,
    directed = FALSE,
    vertices = igraph::as_data_frame(g0, what = "vertices")
  )

  if (remove_isolates) {
    keep2 <- igraph::degree(g) > 0
    g <- igraph::induced_subgraph(g, vids = which(keep2))
  }

  if (igraph::vcount(g) < 2 || igraph::ecount(g) == 0) {
    stop("No drawable graph after filtering/removing isolates.")
  }

  # ---- (D) node frequency from dfm: robust way ----
  freq_all <- quanteda::featfreq(dfm)  # named numeric, safe for dfm S4
  freq <- as.numeric(freq_all[igraph::V(g)$name])
  freq[is.na(freq)] <- 0
  igraph::V(g)$freq <- freq
  igraph::V(g)$freq_plot <- if (isTRUE(use_freq_as_size)) freq else rep(1, length(freq))

  # ---- (E) community detection ----
  if (community == "fast_greedy") {
    com <- igraph::cluster_fast_greedy(g, weights = igraph::E(g)$weight)
    igraph::V(g)$comm <- as.integer(igraph::membership(com))
  } else if (community == "walktrap") {
    com <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight)
    igraph::V(g)$comm <- as.integer(igraph::membership(com))
  } else if (community == "edge_betweenness") {
    com <- igraph::cluster_edge_betweenness(g, weights = igraph::E(g)$weight)
    igraph::V(g)$comm <- as.integer(igraph::membership(com))
  } else {
    igraph::V(g)$comm <- 1L
  }
  igraph::V(g)$comm <- factor(igraph::V(g)$comm)

  # ---- (F) edge linetype: KH Coder rule (FIXED) ----
  # Within-community edges: solid
  # Between-community edges: dotted
  ends <- igraph::ends(g, igraph::E(g), names = TRUE)
  vname <- igraph::V(g)$name
  v_index <- setNames(seq_along(vname), vname)

  a <- v_index[ends[, 1]]
  b <- v_index[ends[, 2]]
  memb <- igraph::V(g)$comm
  same_comm <- memb[a] == memb[b]

  igraph::E(g)$lty <- ifelse(same_comm, "solid", "dotted")

  # ---- (G) tidygraph + plot ----
  tg <- tidygraph::as_tbl_graph(g)

  set.seed(seed)

  size_multiplier <- bubble_size / 100
  node_tbl <- tidygraph::as_tibble(tg, active = "nodes")
  size_aes <- if (all(node_tbl$freq_plot <= 0)) 1 else "freq_plot"

  # Edge layer (linetype identity)
  p <- ggraph::ggraph(tg, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(linetype = lty),
      colour = "grey50",
      linewidth = 0.5
    ) +
    ggplot2::scale_linetype_identity() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_family = font_family)

  # Node layer
  if (!isTRUE(bw_groups)) {
    # Color mode: allow alpha control via fill_alpha
    fill_alpha <- as.numeric(fill_alpha)
    if (is.na(fill_alpha)) fill_alpha <- 0.75
    fill_alpha <- max(min(fill_alpha, 1), 0)

    p <- p +
      ggraph::geom_node_point(
        ggplot2::aes(size = .data[[size_aes]], fill = comm),
        shape = 21,
        colour = "grey30",
        stroke = 0.8,
        alpha = fill_alpha
      )
  } else {
    # Grayscale mode: alpha fixed for readability
    p <- p +
      ggraph::geom_node_point(
        ggplot2::aes(size = .data[[size_aes]], fill = comm),
        shape = 21,
        colour = "grey30",
        stroke = 0.8,
        alpha = 0.95
      )
  }

  # Size scale (+ legend option)
  if (isTRUE(use_freq_as_size)) {
    if (isTRUE(show_size_legend)) {
      p <- p + ggplot2::scale_size_area(
        max_size = 20 * size_multiplier,
        name = "Frequency"
      )
    } else {
      p <- p + ggplot2::scale_size_area(
        max_size = 20 * size_multiplier,
        guide = "none"
      )
    }
  } else {
    p <- p + ggplot2::scale_size_area(
      max_size = 20 * size_multiplier,
      guide = "none"
    )
  }

  # Labels
  if (isTRUE(label_repel)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
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
  } else {
    p <- p + ggraph::geom_node_text(
      ggplot2::aes(label = name),
      family = font_family,
      size = label_size,
      colour = label_colour
    )
  }

  # Fill scale: color vs grayscale
  if (isTRUE(bw_groups)) {
    # Make the darkest shade brighter (avoid near-black)
    p <- p + ggplot2::scale_fill_grey(start = 0.88, end = 0.35, name = "Community")
  } else {
    p <- p + ggplot2::scale_fill_hue(name = "Community")
  }

  info <- list(
    nodes = igraph::vcount(g),
    edges = igraph::ecount(g),
    density = igraph::edge_density(g, loops = FALSE),
    threshold = th_use
  )

  list(graph = g, tbl_graph = tg, plot = p, info = info)
}
