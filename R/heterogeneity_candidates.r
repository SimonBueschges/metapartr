heterogeneity_candidates <-
  function(moderators, model, subset, node.index) {
    Combination.tables <-
      lapply(moderators, function(x) {
        moderator <- model$data[subset, ][[x]]
        unique.moderator <- unique(moderator)
        unique.length.moderator <- length(unique.moderator)

        if (is.factor(moderator) |
            is.character(moderator)) {
          Combinations <- lapply(
            1:unique.length.moderator,
            function(x) {
              combn(unique(moderator), x)
            }
          )

          Combinations <- Combinations[-unique.length.moderator]


          Combination.vectors.list <-
            unlist(lapply(Combinations, function(y) {
              lapply(1:ncol(y), function(z) {
                y[, z]
              })
            }),
            recursive = FALSE
            )


          Table <- data.frame(
            Group.1 = sapply(
              Combination.vectors.list,
              function(x) {
                paste(x, collapse = " + ")
              }
            ),
            Group.2 = sapply(
              sapply(Combination.vectors.list, setdiff, x = unique.moderator),
              function(x) {
                paste(x, collapse = " + ")
              }
            )
          )

          Table[!duplicated(t(apply(Table, 1, sort))), ]
        } else {
          data.frame(Cutpoint = sort(unique(moderator)))
        }
      })

    names(Combination.tables) <- moderators


    # Fit Model for the i - classes
    message("Fitting Models, this takes a moment")


    Models <- lapply(
      moderators,
      function(x) {
        message(paste("Currently for", x))

        if (nrow(Combination.tables[[x]]) > 0) {
          if (all(names(Combination.tables[[x]]) == c("Group.1", "Group.2"))) {
            pb <-
              txtProgressBar(
                style = 3,
                min = 0,
                max = length(Combination.tables[[x]]$Group.1) * 2
              )
            Current.Models <- list(
              Group.1 =
                lapply(Combination.tables[[x]]$Group.1, function(y) {
                  setTxtProgressBar(pb, which(Combination.tables[[x]]$Group.1 == y))
                  update.meta(model, subset = with(model$data, get(x) %in% unlist(strsplit(y, " \\+ "))) &
                                subset)
                }),
              Group.2 =
                lapply(Combination.tables[[x]]$Group.2, function(y) {
                  setTxtProgressBar(
                    pb,
                    which(Combination.tables[[x]]$Group.2 == y) + length(Combination.tables[[x]]$Group.1)
                  )
                  update.meta(model, subset = with(model$data, get(x) %in% unlist(strsplit(y, " \\+ "))) &
                                subset)
                })
            )
            close(pb)
            Current.Models
          } else {
            pb <-
              txtProgressBar(
                style = 3,
                min = 0,
                max = length(Combination.tables[[x]]$Cutpoint[-1]) * 2
              )
            Current.Models <- list(
              Group.1 =
                lapply(Combination.tables[[x]]$Cutpoint, function(y) {
                  if (Combination.tables[[x]]$Cutpoint[1] == y) {
                    NULL
                  } else {
                    setTxtProgressBar(pb, which(Combination.tables[[x]]$Cutpoint[-1] == y))
                    update.meta(model, subset = with(model$data, get(x) < y) &
                                  subset)
                  }
                }),
              Group.2 =
                lapply(Combination.tables[[x]]$Cutpoint, function(y) {
                  if (Combination.tables[[x]]$Cutpoint[1] == y) {
                    NULL
                  } else {
                    setTxtProgressBar(
                      pb,
                      which(Combination.tables[[x]]$Cutpoint[-1] == y) + length(Combination.tables[[x]]$Cutpoint[-1])
                    )
                    update.meta(model, subset = with(model$data, get(x) >= y) &
                                  subset)
                  }
                })
            )
            close(pb)
            Current.Models
          }
        }
      }
    )


    names(Models) <- moderators


    # theta the pooled effect size in the total set of studies.
    # theta_i the pooled effect size in the -ith class
    Combination.tables <-
      lapply(
        moderators,
        function(x) {
          if (nrow(Combination.tables[[x]]) > 0) {
            if (all(names(Combination.tables[[x]]) == c("Group.1", "Group.2"))) {
              Combination.tables[[x]]$Q.B <-
                sapply(
                  1:length(Models[[x]][[1]]),
                  function(i) {
                    sum(
                      ((Models[[x]][[1]][[i]]$TE.fixed - model$TE.fixed)^2 * 1 / Models[[x]][[1]][[i]]$seTE.fixed^
                         2
                      ),
                      ((Models[[x]][[2]][[i]]$TE.fixed - model$TE.fixed)^2 * 1 / Models[[x]][[2]][[i]]$seTE.fixed^
                         2
                      )
                    )
                  }
                )
            } else {
              Combination.tables[[x]]$Group.1 <-
                paste(" < ", Combination.tables[[x]]$Cutpoint)
              Combination.tables[[x]]$Group.2 <-
                paste(" >= ", Combination.tables[[x]]$Cutpoint)
              Combination.tables[[x]]$Cutpoint <- NULL

              Combination.tables[[x]]$Q.B <-

                  sapply(
                    1:length(Models[[x]][[1]]),
                    function(i) {
                      sum(
                        ((Models[[x]][[1]][[i]]$TE.fixed - model$TE.fixed)^2 * 1 / Models[[x]][[1]][[i]]$seTE.fixed^
                           2
                        ),
                        ((Models[[x]][[2]][[i]]$TE.fixed - model$TE.fixed)^2 * 1 / Models[[x]][[2]][[i]]$seTE.fixed^
                           2
                        )
                      )
                    }
                  )

            }

            # theta_ij is the effect size in the -jth study in the ith_class
            Combination.tables[[x]]$Q.W <-
              Model$Q - Combination.tables[[x]]$Q.B

            Combination.tables[[x]]$pval.Q.B <-
              pchisq(Combination.tables[[x]]$Q.B, (2 - 1), lower.tail = FALSE)

            Combination.tables[[x]]$logworth <-
              -log10(Combination.tables[[x]]$pval.Q.B)

            # Ratio
            # 2 being the number of splits
            Combination.tables[[x]]$Ratio.QB.QW <-
              (Combination.tables[[x]]$Q.B / (2 - 1)) /
              (Combination.tables[[x]]$Q.W / (model$k - 2))

            Combination.tables[[x]]
          }
        }
      )

    names(Combination.tables) <- moderators

    list(Models = Models, Combinations.tables = Combination.tables)
  }