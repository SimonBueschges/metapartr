heterogeneity_candidates <- function(moderators, model, unique.length.moderators) {

  browser()

  Combinations.list <-
    lapply(moderators, function(x) {

      moderator <- model$data[[x]]
      unique.length.moderator <- unique.length.moderators[[x]]

      if (is.factor(moderator) |
          is.character(moderator)) {

        Combinations <- lapply(
          1:unique.length.moderator,
          function(x) {
            combn(unique(moderator), x)
          }
        )

        Combinations[-unique.length.moderator]
      }
    })

  names(Combinations.list) <- moderators

  Combinations.list <-
    lapply(Combinations.list, function(x) {
      lapply(x, function(y) {
        lapply(1:ncol(y), function(z) {
          y[, z]
        })
      })
    })


  Q_H <- model$Q

  # Fit Model for the i - classes

  I.classes <- lapply(Combinations.list, unlist, recursive = FALSE)

  Models <- lapply(
    names(I.classes),
    function(x) {
      lapply(I.classes[[x]],
             function(y){
               update.meta(model, subset = with(model$data, get(x) == y))
             }
      )
    }
  )

  names(Models) <- moderators

  Models <-  Models[sapply(Models,length) > 0]
  # theta the pooled effect size in the total set of studies.
  # theta_i the pooled effect size in the -ith class
  Q_B <-
    lapply(
      1:length(Models),
      function(x){
        sum(
          sapply(
            1:length(x),
            function(i) {
              (Models[[x]][[i]]$TE.fixed - model$TE.fixed)^2 * 1 / Models[[x]][[i]]$seTE.fixed^2
            }
          )
        )
      }
    )

  # theta_ij is the effect size in the -jth study in the ith_class
  Q_W <- Model$Q - Q_B

  # Ratio
  # 2 being the number of splits
  (Q_B / (2 - 1)) / (Q_W / nrow(Model$data) - 2)


  # meta:::hetcalc(
  #   model$TE, seTE = model$seTE,
  #   method.tau = model$method.tau,
  #   TE.tau = model$TE.tau,
  #   method.tau.ci =  model$method.tau.ci, level.hetstats = 0.95)$Q

  # data.frame(
  #   Moderator,
  #   Q_B,
  #   Logworth,
  #   Q_W
  #   Ratio.of.QB.and.QW.
  #   Cutpoint
  # )

}