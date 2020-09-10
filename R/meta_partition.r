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
#' An object of class \code{c("metagen", "meta")} with corresponding
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
meta_partition <- function(model,
                           moderators,
                           n.to.small = 8,
                           c = 1,
                           auto = TRUE,
                           control = rpart.control(
                             xval = 10,
                             minbucket = 3,
                             minsplit = 6,
                             cp = 0.0001
                           )) {
  # Check arguments ---------------------------------------------------------

  # Needs to check model class to be "meta"
  assert_class(model, classes = "meta")

  # Needs to check that that is part of model
  assert_data_frame(model$data)

  # Needs to check that the moderators are part of the model data

  assert_subset(moderators, names(model$data))


  assert_logical(auto)


  # Variable Definitions ----------------------------------------------------



  # create formula for rpart
  formula <-
    as.formula(paste(
      "TE ~",
      paste(paste0(
        "`", moderators, "`"
      ), collapse = "+")
    ))

  # define unique length of moderators
  unique.length.moderators <-
    vapply(moderators,
      FUN.VALUE = 1L,
      function(x) {
        length(unique(model$data[[x]]))
      }
    )

  if (any(unique.length.moderators <= 1)) {
    warning("moderators are only reasonable with more than one unique value")
  }

  Nodes <- list(list(Parent = NA, Terminal.Node = FALSE, Model = model))
  node.index <- 0

  # Recursion ----------------------------------------------------------------

  # The approach consists in assessing heterogeneity
  # of a sample of effect sizes and partitioning this variability
  # with an algorithm that finds the best moderator of each partition.



  # > Step 1: study of heterogeneity ----------------------------------------


  # step 1 is to test the initial hypothesis of
  # one real effect size with the fixed effect
  # model, for testing homogeneity.

  # Q_H is used to test homogeneity
  # sum((Model$TE-Model$TE.fixed)^2 * 1/Model$seTE^2)


  while (node.index < length(Nodes)) {
    # If the set of effect sizes showed a significant heterogeneity, we
    # then calculated the amount of heterogeneity with I2
    message(
      paste(
        "Partitioning\n",
        "Q.H =",
        round(Nodes[[node.index + 1]][["Model"]]$Q, 2),
        ",",
        "pval.Q.H =",
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
        Nodes[[node.index + 1]][["Model"]]$k > n.to.small) {
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


    # > Step 2: partition analysis  ------------------------------------------



    # Afterwards, we will try to explain
    # this heterogeneity by the moderators
    # that were established a priori in the design of meta-analysis.

    # Second step consists of dividing the subset of effect sizes by the
    # moderator that explains the largest amount
    # of variability between resulting subsets in relation
    # with variability within each of the potential subsets

    # Thus, a partition method is proposed to explore the importance
    # of each moderator on explaining heterogeneity. This partition
    # is based on classification and regression trees

    # The method starts from a
    # set of M effect sizes and a set of moderators (X1, X2, . . ., Xp)

    # M <- length(Nodes[[node.index +1]][["Model"]]$TE)
    # p <- length(moderators)

    # The aim of the algorithm of partition is to find r disjoint
    # classes (m1, m2. . .mr) from the set of effect sizes such that the
    # classes will be the most homogeneous within them,
    # while being the most heterogeneous
    # between them.


    # The algorithm of this procedure is:
    #
    # (1) Looking for the lowest number of classes for each
    # moderator that maximizes Q_B value (usually two).
    # The combination of levels of the moderator will depend on the type of moderator:
    #
    # if it is nominal, all possible ways to combine the categories would
    # need to be evaluated;


    # if it is ordinal, the combination of categories must preserve the
    # order in the data;
    # if the moderator is continuous, the algorithm finds the value that maximizes
    # the interclass distance, and splits the moderator into two groups
    # regarding this value.

    # (2) The first partition is chosen to maximize the difference in
    # the response (effect size) between the levels of the moderator



    # The technique involves the partition of the sum of squares from the test of
    # homogeneity into two components Q_H = Q_B + Q_W



    # In each partition, this method will minimize the intraclass distance (QW) and so maximize the
    # interclass distance (QB). An indicator of importance of the each factor is the ratio of both
    # distances.

    if (keep.partitioning) {
      if (is.null(Nodes[[node.index + 1]][["Model"]]$subset)) {
        subset <- rep(TRUE, Nodes[[node.index + 1]][["Model"]]$k)
      } else {
        subset <- Nodes[[node.index + 1]][["Model"]]$subset
      }

      moderator.length.more.than.one <-
        vapply(
          moderators, function(x) {
            moderator <- Nodes[[node.index + 1]][["Model"]]$data[subset, ][[x]]
            unique.moderator <- unique(moderator)
            unique.length.moderator <- length(unique.moderator)
            unique.length.moderator > 1
          },
          logical(1)
        )

      if (any(moderator.length.more.than.one)) {
        Candidates <-
          heterogeneity_candidates(
            moderators = moderators,
            model = Nodes[[node.index + 1]][["Model"]],
            subset = subset,
            node.index =node.index
            )

        Top.Candidates <- do.call(
          "rbind",
          lapply(
            Candidates[["Combinations.tables"]],
            function(x) {
              Candidate.Combination.entry <- x[which.max(x$Q.B), ]
              cbind(Split.ID = which.max(x$Q.B), Candidate.Combination.entry)
            }
          )
        )

        message("Top Results per Moderator:")

        print(Top.Candidates[order(Top.Candidates$Q.B), ], digits = 4)

        if (isTRUE(auto) | isFALSE(auto)) {
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
        # message(paste("Press enter to choose",
        #               row.names(Top.Candidates[which.max(Top.Candidates$Q.B),]),
        #               "for the split,\n or enter the name of a moderator to view their detailed results"
        # )
        # )
        # split.moderator <- NA
        # while (is.na(split.moderator)) {
        #   answer <- readline("Keep Partitioning? y/n: ")
        #   keep.partitioning <-
        #     if (answer %in% c("y", "yes", "Y", "YES", "Yes", "TRUE", "T", 1)) {
        #       TRUE
        #     } else if (answer %in% c("n", "no", "N", "NO", "No", "FALSE", "F", 0)) {
        #       FALSE
        #     } else {
        #       NA
        #     }
        # }
      }


      Nodes <-
        append(
          Nodes,
          list(list(
            Parent = node.index + 1,
            Model = Candidates[["Models"]][[selected.moderator]][["Group.1"]][[selected.split]]
          ))
        )

      Nodes <-
        append(
          Nodes,
          list(list(
            Parent = node.index + 1,
            Model = Candidates[["Models"]][[selected.moderator]][["Group.2"]][[selected.split]]
          ))
        )
      Nodes[[node.index + 1]][["Terminal.Node"]] <- FALSE
      } else {
        Nodes[[node.index + 1]][["Model"]] <- TRUE
      }
    } else {
      Nodes[[node.index + 1]][["Terminal.Node"]] <- TRUE
    }

    # (3) Testing if the classes obtained in the first partition are
    # homogeneous using QH and I2.

    # (4) If any class is heterogeneous,
    # then the step 2 will be repeated with the rest of moderators



    # Then, step 2 will be repeated until we reach one of these three situations:
    # (1) to reach a subset of effect sizes that
    # is homogeneous under the fixed effect model,
    # (2) to reach a subset of effect sizes that is still heterogeneous
    # but that it is too small to keep partitioning, or
    # (3) to reach a subset of effect sizes that is still heterogeneous
    # but the considered moderators are not able to
    # explain the heterogeneity in this subset.


    # (5) The partition process will stop when all classes are homogeneous,
    # when they are final due to small sample size,
    # or if none of the moderator explains the heterogeneity

    node.index <- node.index + 1

    message(paste("Finished partitioning node", node.index, "\n\n"))
  } # End while



  # Step 3: integrating effect sizes. ---------------------------------------


  Nodes
}
