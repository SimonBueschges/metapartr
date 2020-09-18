# Input meta object (metagen()) and a priori defined
# moderators (vector of variable names)
# Partition algorithm and options
# There are two criteria to do the partition: the p-values of QB or the
# ratio (QB/(p-1))/(QW/(k-p)).
#
# minimum number of studies to stop partitioning n.min = 8



#' Meta Partition
#'
#' Does Meta-Partition as discribed in ..
#'
#' @param model sdfs
#' @param moderators sfssdf
#' @param c sdfdf
#' @param control sdf
#' @param ...
#'
#' @return
#' An object of class \code{"meta"} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#'
#' \item{tau2}{Between-study variance \eqn{\tau^2}.}
#' \item{se.tau2}{Standard error of \eqn{\tau^2}.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{zval.random.w, pval.random.w}{z-value or t-value and
#'   corresponding p-value for test of treatment effect in subgroups
#'   (random effects model) - if \code{byvar} is not missing.}
#'
#'
#' @examples
#' data(Method.Study.Data)
#' Model <-
#'   meta::metacor(
#'     ESr, n,
#'     data = Method.Study.Data,
#'     studlab = Study,
#'     comb.random = FALSE
#'   )
#' @author Simon Büschges
#'
#' @references Ortega Z, Martín-Vallejo J, Mencía A,
#' Galindo-Villardón MP, Pérez-Mellado V (2016)
#' Introducing Meta-Partition, a Useful Methodology to
#' Explore Factors That Influence Ecological Effect
#' Sizes. PLoS ONE 11(7): e0158624. doi:10.1371/
#'   journal.pone.0158624
#'
#'
#' @export
#'
#' @import stringr data.table checkmate
meta_partition <- function(model,
                           moderators,
                           n.too.small = 8,
                           auto = TRUE) {
  # Check arguments ---------------------------------------------------------

  # Check model class contains "meta"
  assert_class(model, classes = "meta")
  # Needs to check that that is part of model
  assert_data_frame(model$data)
  # Needs to check that the moderators are part of the model data
  assert_subset(moderators, names(model$data))

  assert_logical(auto)


  unique.length.moderators <-
    vapply(moderators,
           FUN.VALUE = 1L,
           function(x) {
             length(unique(model$data[[x]]))
           })

  if (any(unique.length.moderators <= 1)) {
    warning("moderators are only reasonable with more than one unique value")
  }

  Nodes <- list(
    list(
      Parent = NA,
      Terminal.Node = FALSE,
      Model = model,
      Combinations.tables = NA,
      choosen.split = NA,
      splits = ""
    )
  )

  groups <- rep(NA, nrow(model$data))

  node.index <- 0


  while (node.index < length(Nodes)) {
    # If the set of effect sizes showed a significant heterogeneity, we
    # then calculated the amount of heterogeneity with I2
    message(
      paste(
        "Evaluting Heterogeneity\n",
        "Q.T =",
        round(Nodes[[node.index + 1]][["Model"]]$Q, 2),
        ",",
        "pval.Q.T =",
        round(Nodes[[node.index + 1]][["Model"]]$pval.Q, 5),
        "\n",
        "I^2 =",
        round(Nodes[[node.index + 1]][["Model"]]$I2, 2),
        "\n",
        "Sample.Size =",
        Nodes[[node.index + 1]][["Model"]]$k
      )
    )

    if (auto == TRUE) {
      if (Nodes[[node.index + 1]][["Model"]]$pval.Q < 0.1 &
          Nodes[[node.index + 1]][["Model"]]$k > n.too.small) {
        keep.partitioning <- TRUE
      } else {
        keep.partitioning <- FALSE
      }
    } else {
      keep.partitioning <- NA
      while (is.na(keep.partitioning)) {
        answer <- readline("Keep Partitioning? y/n: ")
        keep.partitioning <-
          if (answer %in% c("y", "yes", "Y", "YES", "Yes", "TRUE", "T", 1)) {
            TRUE
          } else if (answer %in% c("n", "no", "N", "NO", "No", "FALSE", "F", 0)) {
            FALSE
          } else {
            NA
          }
      }
    }

    if (keep.partitioning) {
      if (is.null(Nodes[[node.index + 1]][["Model"]]$subset)) {
        subset <- rep(TRUE, Nodes[[node.index + 1]][["Model"]]$k)
      } else {
        subset <- Nodes[[node.index + 1]][["Model"]]$subset
      }

      moderator.length.more.than.one <-
        vapply(moderators, function(x) {
          moderator <- Nodes[[node.index + 1]][["Model"]]$data[subset, ][[x]]
          unique.moderator <- unique(moderator)
          unique.length.moderator <- length(unique.moderator)
          unique.length.moderator > 1
        },
        logical(1))

      if (any(moderator.length.more.than.one)) {
        Candidates <-
          heterogeneity_candidates(
            moderators = moderators,
            model = Nodes[[node.index + 1]][["Model"]],
            subset = subset,
            node.index = node.index
          )

        Top.Candidates <- do.call("rbind",
                                  lapply(Candidates[["Combinations.tables"]],
                                         function(x) {
                                           Candidate.Combination.entry <- x[which.max(x$Q.B), ]
                                           cbind(Split.ID = which.max(x$Q.B), Candidate.Combination.entry)
                                         }))

        message("Top Results per Moderator:")

        print(Top.Candidates[order(Top.Candidates$Q.B), ], digits = 4)

        if (isTRUE(auto)) {
          selected.moderator <-
            row.names(Top.Candidates[which.max(Top.Candidates$Q.B), ])
          selected.split <-
            Top.Candidates[selected.moderator, "Split.ID"]


          message(
            paste0(
              "Using the ",
              selected.split,
              "-th possible split of ",
              selected.moderator,
              " for the split"
            )
          )
        } else {
          repeat.message <- TRUE

          while (isTRUE(repeat.message)) {
            message(
              "Select the moderator and Split.ID for the partition,
            by typing the name followed by the Split.ID. To view the details of
            a moderator type the name of the moderator followed
                    by the word details:"
            )

            answer.selection <- readline("")

            if (!str_detect(answer.selection,
                            pattern = paste0(rownames(Top.Candidates), collapse = "|"))) {
              message("Please type a valid moderator name")
            } else{
              if (!str_detect(answer.selection,
                              pattern = paste0(c(
                                "details",
                                unique(unlist(
                                  lapply(Candidates$Combinations.tables,
                                         function(x) {
                                           1:nrow(x)
                                         }),
                                  use.names = FALSE
                                ))
                              ),
                              collapse = "|"))) {
                message("Please type details or a valid Split.ID")
              } else {
                if(str_detect(answer.selection,
                              pattern = "details"
                )){
                  lapply(names(Candidates$Combinations.tables),
                         function(x){
                           if(str_detect(answer.selection,
                                         pattern = x
                           )){
                             print(Candidates$Combinations.tables[[x]])
                           }
                         }
                         )
                } else {

                  selected.moderator <-
                    str_extract(answer.selection, paste0(rownames(Top.Candidates), collapse = "|"))

                  if (!str_detect(answer.selection,
                                  pattern = paste0(
                                    1:nrow(Candidates$Combinations.tables[[selected.moderator]]),
                                  collapse = "|"))) {

                    message("Please type a valid Split.ID")
                  } else {

                  selected.split <-
                    as.numeric(
                      str_extract(answer.selection, "[[:digit:]]+")
                    )

                  repeat.message <- FALSE


                  message(
                    paste0(
                      "Using the ",
                      selected.split,
                      "-th possible split of ",
                      selected.moderator,
                      " for the split"
                    )
                  )
                  }
                }
              }
            }
          }
        }

          Nodes <-
            append(Nodes,
                   list(
                     list(
                       Parent = node.index + 1,
                       Model = Candidates[["Models"]][[selected.moderator]][["Group.1"]][[selected.split]],
                       splits = paste(Nodes[[node.index + 1]][["splits"]],
                                      "->",
                                      selected.moderator,
                                      Candidates[["Combinations.tables"]][[selected.moderator]][selected.split, ][["Group.1"]])
                     )
                   ))

          Nodes <-
            append(Nodes,
                   list(
                     list(
                       Parent = node.index + 1,
                       Model = Candidates[["Models"]][[selected.moderator]][["Group.2"]][[selected.split]],
                       splits = paste(Nodes[[node.index + 1]][["splits"]],
                                      "->",
                                      selected.moderator,
                                      Candidates[["Combinations.tables"]][[selected.moderator]][selected.split, ][["Group.2"]])
                     )
                   ))

          Nodes[[node.index + 1]][["Terminal.Node"]] <- FALSE

          Nodes[[node.index + 1]][["Combinations.tables"]] <-
            Candidates$Combinations.tables

          Nodes[[node.index + 1]][["choosen.split"]] <-
            Candidates[["Combinations.tables"]][[selected.moderator]][selected.split, ]
        } else {
          Nodes[[node.index + 1]][["Model"]] <- TRUE
        }
      } else {
        Nodes[[node.index + 1]][["Terminal.Node"]] <- TRUE
        groups[Nodes[[node.index + 1]]$Model$subset] <-
          node.index + 1
      }


      node.index <- node.index + 1

      message(paste("Finished partitioning node", node.index, "\n\n"))
    } # End while

    list(Nodes = Nodes, groups = groups)
  }
