

#' Title
#'
#' @param DE_results_comb
#'
#' @returns
#' @export
#'
#' @examples
plot_L2FC_scatter <- function(
    DE_Results_comb,
    x=NULL,
    y=NULL,
    genes="All",
    flavor="Min",
    fdr=0.05,
    l2fc=0,
    symmetric=TRUE,
    guides=TRUE,
    pt.size=1,
    annotate_genes=NULL
    ){

  axes = select_comb_DE_results(DE_Results_comb, X=x, Y=y, type = "L2FC")
  lm = max(
    abs(DE_Results_comb[[axes[1]]]),
    abs(DE_Results_comb[[axes[2]]]),
    na.rm = TRUE
  )
  p1 = ggplot2::ggplot() +

    {
      if(guides){
        ggplot2::geom_hline(yintercept = 0, color="red")

      }else{
        NULL
      }
    } +
    {
      if(guides){
        ggplot2::geom_vline(xintercept = 0, color="red")
      }else{
        NULL
      }
    } +

    ggplot2::aes(
      x=.data[[axes[1]]],
      y=.data[[axes[2]]],
      label = ifelse(Gene %in% annotate_genes, Gene, NA)
    ) +

    {if(genes %in% c("All", "X", "Only_one")) {
      ggplot2::geom_point(
      data = DE_Results_comb[
        which(
          DE_Results_comb[[stringr::str_replace(axes[1], "L2FC", "QVAL")]] >= fdr |
            (DE_Results_comb[[stringr::str_replace(axes[1], "L2FC", "QVAL")]] < fdr &
               abs(DE_Results_comb[[axes[1]]]) <= l2fc)
        ),
      ],
      pch=switch(
        flavor,
        Min = ".",
        Classic=16,
        Fill=21,
        21
      ),
      size=pt.size,
      fill="blue"
    )}else{
      NULL
      }
    } +

    {
      if(!is.null(annotate_genes)){
        ggrepel::geom_text_repel(data = DE_Results_comb[
          which(
            DE_Results_comb[[stringr::str_replace(axes[1], "L2FC", "QVAL")]] >= fdr |
              (DE_Results_comb[[stringr::str_replace(axes[1], "L2FC", "QVAL")]] < fdr &
                 abs(DE_Results_comb[[axes[1]]]) <= l2fc)
          ),
        ])
      }else{
        NULL
      }
    } +

    {if(genes %in% c("All", "Y", "Only_one")){
      ggplot2::geom_point(
        ggplot2::aes(
          x=.data[[axes[1]]],
          y=.data[[axes[2]]]
        ),
        data = DE_Results_comb[
          which(
            DE_Results_comb[[stringr::str_replace(axes[2], "L2FC", "QVAL")]] >= fdr |
              (DE_Results_comb[[stringr::str_replace(axes[2], "L2FC", "QVAL")]] < fdr &
                 abs(DE_Results_comb[[axes[2]]]) <= l2fc)
          ),
        ],
        pch=switch(
          flavor,
          Min = ".",
          Classic=16,
          Fill=21,
          21
        ),
        size=pt.size,
        fill="green"
      )
    }else{
      NULL
      }
    } +

    {
      if(!is.null(annotate_genes)){
        ggrepel::geom_text_repel(data = DE_Results_comb[
          which(
            DE_Results_comb[[stringr::str_replace(axes[2], "L2FC", "QVAL")]] >= fdr |
              (DE_Results_comb[[stringr::str_replace(axes[2], "L2FC", "QVAL")]] < fdr &
                 abs(DE_Results_comb[[axes[2]]]) <= l2fc)
          ),
        ])
      }else{
        NULL
      }
    } +

    {
      if(genes %in% c("All", "Common", "X", "Y")){
        ggplot2::geom_point(
          ggplot2::aes(
            x=.data[[axes[1]]],
            y=.data[[axes[2]]]
          ),
          data = DE_Results_comb[
            which(
              (DE_Results_comb[[stringr::str_replace(axes[1], "L2FC", "QVAL")]] < fdr &
                 abs(DE_Results_comb[[axes[1]]]) > l2fc) &
                (DE_Results_comb[[stringr::str_replace(axes[2], "L2FC", "QVAL")]] < fdr &
                   abs(DE_Results_comb[[axes[2]]]) > l2fc)

            ),
          ],
          pch=switch(
            flavor,
            Min = ".",
            Classic=16,
            Fill=21,
            21
          ),
          size=pt.size,
          fill="red"
        )
      }else{
        NULL
      }
    } +

    {
      if(!is.null(annotate_genes)){
        ggrepel::geom_text_repel(
          data = DE_Results_comb[
            which(
              (DE_Results_comb[[stringr::str_replace(axes[1], "L2FC", "QVAL")]] < fdr &
                 abs(DE_Results_comb[[axes[1]]]) > l2fc) &
                (DE_Results_comb[[stringr::str_replace(axes[2], "L2FC", "QVAL")]] < fdr &
                   abs(DE_Results_comb[[axes[2]]]) > l2fc)

            ),
          ]
        )
      }else{
        NULL
      }
    } +

    ggplot2::scale_x_continuous(limits = if(symmetric) c(-lm, lm) else NULL) +
    ggplot2::scale_y_continuous(limits = if(symmetric) c(-lm, lm) else NULL) +

    xlab(
      paste0(
        "\n",
        stringr::str_replace(axes[1], "_", " ")

      )
    ) +

    ylab(
      paste0(
        stringr::str_replace(axes[2], "_", " "),
        "\n"
      )
    )


  return(p1)



}
