
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(patchwork))

qc_plot <- function(data, group, orders, filename,n_gene=50,addEllipses=FALSE, pca_legend="Group", ellipse_size = 0.8) {
    # 筛选高变SD基因
    hvg <- function(data, n_gene) {
        tops <- apply(data, 1, sd) %>%
            sort(decreasing = T) %>%
            .[1:n_gene] %>%
            names()
        return(data[tops, ])
    }

    # PCA plot
    pca <- FactoMineR::PCA(t(data), scale.unit = TRUE, ncp = 2, graph = F)

    p1 <- factoextra::fviz_pca_ind(
        X = pca, ## pca对象
        axes = 1:2, ## 展示的两个主成分
        geom = "point", ## 展示individual的形式
        habillage = factor(group), ## individual用来分组的变量
        legend.title = pca_legend, ## 分组变量的title
        palette = "set2", ## 颜色面板
        addEllipses = addEllipses, ## 是否绘制椭圆
        ellipse.level = ellipse_size, ## 椭圆的大小
        title = "PCA plot", ## 标题
        mean.point = FALSE ## 删除每个组的重心
    ) +
        theme(plot.title = element_text(hjust = 0.5, size = 20))

    # 第二张图：聚类树图

    p2 <- hvg(data, 5000) %>%
        t() %>%
        dist() %>%
        hclust() %>%
        ggdendrogram()

    # 第三张图：聚类热图
    mat <- hvg(data,n_gene) %>% 
        t() %>% 
        scale() %>% 
        t()
    p3 <- Heatmap(mat,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_order = orders, # 样品顺序
        row_order = NULL,
        column_title = "",
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        show_heatmap_legend = TRUE,
        show_row_names=ifelse(n_gene>50,FALSE,TRUE),
        heatmap_legend_param = list(title = "")
    )
    # 拼图
    p123 <- (p1 + p2) / wrap_elements(plot = grid.grabExpr(draw(p3)), clip = T) +
        plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A")
    pdf(file=filename,width=12, height=12)
    print(p123)
    dev.off()
}

dt <- data.table::fread(snakemake@input[[1]], data.table = FALSE) %>% 
    dplyr::select(Symbol,starts_with("tpm"))

unique_exprs <- function(exprs_df,  method = "median") {
    exprs_df %>%
    dplyr::mutate(ref = apply(across(where(is.numeric)), 1, method)) %>%
    dplyr::arrange(desc(ref)) %>% # 把表达量的平均值按从大到小排序
    dplyr::select(-ref) %>% # 反向选择去除rowMean这一列
    dplyr::distinct(Symbol, .keep_all = T)
}

dt %<>% unique_exprs %>% 
    tibble::column_to_rownames("Symbol") 

colnames(dt) %<>% stringr::str_remove("tpm_")

sample_order <- snakemake@params[["sample_order"]]
group <- snakemake@params[["group"]]

if(all(snakemake@params[["group"]]!="")){
    pca_legend <- "Group"
    sample_order <- tidyr::tibble(Sample=sample_order,Group=group) %>% 
        dplyr::arrange(Group,Sample) %>%
        dplyr::pull(Sample)
} else {
    group <- sample_order
    pca_legend <- "Sample"
}

qc_plot(data=log2(dt+1),
    group = group,
    n_gene=snakemake@params[["n_gene"]],
    pca_legend = pca_legend,
    filename = snakemake@output[[1]],
    orders = sample_order,
    addEllipses = snakemake@params[["ellipses"]],
    ellipse_size = snakemake@params[["ellipse_size"]]
)

