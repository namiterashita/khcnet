library(khcnet)
library(quanteda)
library(tidyverse)

corp <- corpus(c(
  "男女平等とフェミニズム",
  "反フェミニズム運動",
  "ジェンダー平等政策"
))



dfm <- dfm(tokens(corp, remove_punct = TRUE))

res_col <- khcnet(
  dfm,
  bw_groups = FALSE,
  fill_alpha = 0.55,
  font_family = "Noto Sans CJK JP",
  show_size_legend = TRUE,
  show_linetype_legend = FALSE,
  bubble_size = 100,
  label_colour = "black",
  label_repel = TRUE
)

# ---- 実行例C：白黒（コミュニティ間を実線にする）
res_bw_solid_inter <- khcnet(
  dfm,
  bw_groups = TRUE,
  bw_intercommunity_linetype = "solid",
  font_family = "sans",
  label_colour = "black",
  label_repel = TRUE
)

res_bw <- khcnet(
  dfm,
  bw_groups = TRUE,
  bw_intercommunity_linetype = "dotted",
  font_family = "sans",
  label_colour = "black",
  label_repel = TRUE
)


res_bw$plot+
  ggplot2::guides(fill = ggplot2::guide_legend(title = NULL))

print(res_bw_solid_inter$plot)

print(res_col$plot)

res_col$plot  +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )


res <- khcnet(
  dfm,
  community = "fast_greedy",
  edges = 60,
  bw_groups = FALSE,
  font_family = "Noto Sans CJK JP",
  label_colour = "gray20"
)
res$plot

res <- khcnet(
  dfm,
  bw_groups = TRUE,
  bw_intercommunity_linetype = "solid"
)
res$plot
table(igraph::E(res$graph)$lty)


res1 <- khcnet(dfm, bw_groups = TRUE, bw_intercommunity_linetype = "dotted")
res2 <- khcnet(dfm, bw_groups = TRUE, bw_intercommunity_linetype = "solid")

table(igraph::E(res1$graph)$lty)
table(igraph::E(res2$graph)$lty)
res <- khcnet(dfm, bw_groups = TRUE, bw_intercommunity_linetype = "solid")

# tidygraph 側の edge データに lty が入ってるか
ed <- tidygraph::as_tibble(res$tbl_graph, active = "edges")
names(ed)
table(ed$lty)
