

#' Title
#'
#' @param DE_results_comb
#'
#' @returns
#' @export
#'
#' @examples
plot_L2FC_scatter <- function(DE_Results_comb, x=NULL, y=NULL, genes="All", flavor="Min", fdr=0.05, l2fc=0){

  axes = select_comb_DE_results(DE_Results_comb, X=x, Y=y, type = "L2FC")

  p1 = ggplot2::ggplot() +

    {if(genes %in% c("All", "X", "Only_one")) {
      ggplot2::geom_point(
      ggplot2::aes(
        x=.data[[axes[1]]],
        y=.data[[axes[2]]]
      ),
      data = DE_Results_comb[
        which(
          DE_Results_comb[[stringr::str_replace(axes[1], "L2FC", "QVAL")]] >= fdr |
            (DE_Results_comb[[stringr::str_replace(axes[1], "L2FC", "QVAL")]] < fdr &
               abs(DE_Results_comb[[axes[1]]]) <= l2fc)
        ),
      ],
      pch=21,
      fill="blue"
    )}else{
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
        pch=21,
        fill="green"
      )
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
          pch=21,
          fill="red"
        )
      }else{
        NULL
      }
    }


  return(p1)



}
