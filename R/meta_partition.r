

# Input meta object (metagen()) and a priori defined
# moderators (vector of variable names)
# Partition algorithm and options
# There are two criteria to do the partition: the p-values of QB or the
# ratio (QB/(p-1))/(QW/(k-p)).
#
# minimum number of studies to stop partitioning n.min = 8

# Step 1

# Test for homogeneity

# Calculate Heterogeneity


# Step 2: partition analysis

# Second step consists of dividing the subset of effect sizes by the
# moderator that explains the largest amount of variability between resulting subsets in relation
# with variability within each of the potential subsets

# classes will be the most homogeneous within them, while being the most heterogeneous
# between them. The technique involves the partition of the sum of squares from the test of
# homogeneity into two components


# QH = QB + QW
# in where:
#   QB = sum_from_i=1_to_r of (theta_hat_i  - theta_hat)^2 * w_i
# and
#   QW = sum_from_i=1_to_r of  sum_from_j=1_to_k (theta_hat_ij  - theta_hat)^2 * w_ij

# An indicator of importance of the each factor is the ratio of both
# distances.


#' Title
#'
#' @param model
#' @param moderators
#' @param c
#' @param control
#' @param ...
#'
#' @return
#' @export
#'
#' @examples data(Method.Study.Data)
#'Model <-
#'meta::metacor(
#'ESr, n,
#' data = Method.Study.Data,
#' studlab = Study,
#' comb.random = FALSE
#' )
Metapart <-
  meta_partition(
    model = Model,
    moderators = c(
      "Study.Type", "Sampling.Effort", "Patch.Area",
      "Home.Range.(ha)", "Length.(cm)",  "Repro","Taxa"
    ),
    c = 1,
    control = rpart.control(
      xval = 10,
      minbucket = 3,
      minsplit = 6,
      cp = 0.0001
    )
  )
meta_partition <- function(model,
                           moderators,
                           n.to.small = 8,
                           c = 1,
                           control = rpart.control(
                             xval = 10,
                             minbucket = 3,
                             minsplit = 6,
                             cp = 0.0001
                           )) {
  browser()
  # Check arguments ---------------------------------------------------------

  # Needs to check model class to be "meta"
  assert_class(model, classes = "meta")

  formula <- as.formula(
    paste(
      "TE ~",
      paste(
        paste0("`",moderators, "`"), collapse = "+"
      )
    )
  )


  # Recursion ----------------------------------------------------------------

  # The approach consists in assessing heterogeneity
  # of a sample of effect sizes and partitioning this variability
  # with an algorithm that finds the best moderator of each partition.



  # > Step 1: study of heterogeneity ----------------------------------------


  # step 1 is to test the initial hypothesis of
  # one real effect size with the fixed effect
  # model, for testing homogeneity.

  # Q_H is used to test homogeneity

  model$Q # sum((Model$TE-Model$TE.fixed)^2 * 1/Model$seTE^2)

  # while
  if(model$pval.Q < 0.1) {


    # If the set of effect sizes showed a significant heterogeneity, we
    # then calculated the amount of heterogeneity with I2

    print(paste("Q_H =", model$Q, "I^2 =", model$I2))





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

    model$TE
    moderators

    M <- length(model$TE)
    p <- length(moderators)

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


    # Data.tree <-
    #   cbind(model$data, TE =  model$TE, wts =  1/model$seTE)
    #
    # tree <- rpart(
    #   formula,
    #   weights = wts,
    #   data = Data.tree,
    #   control = rpart.control(
    #     xval = 10,
    #     minbucket = 3,
    #     minsplit = 6,
    #     cp = 0.0001,
    #     maxdepth = 1
    #   )
    # )

    # The technique involves the partition of the sum of squares from the test of
    # homogeneity into two components Q_H = Q_B + Q_W

    Q_H <- model$Q

    # Fit Model for the i - classes


    split.var <- "Taxa"
    i.classes <- list(
      c("a", "r"),
      c("m", "b")
    )

    Models <- lapply(
      i.classes,
      function(x) {
        update.meta(Model, subset = with(Model$data, Taxa == "a"))
      }
    )
    # theta the pooled effect size in the total set of studies.
    # theta_i the pooled effect size in the -ith class
    Q_B = sum(
      sapply(1:i,
             function(i){
               (Models[[i]]$TE.fixed-model$TE.fixed)^2 * 1/Models[[i]]$seTE.fixed^2
             }
      )
    )

    # theta_ij is the effect size in the -jth study in the ith_class
    Q_W = Model$Q - Q_B

    # Ratio
    # 2 being the number of splits
    (Q_B/(2-1))/(Q_W/nrow(Model$data)- 2)

    # In each partition, this method will minimize the intraclass distance (QW) and so maximize the
    # interclass distance (QB). An indicator of importance of the each factor is the ratio of both
    # distances.



    # (3) Testing if the classes obtained in the first partition are
    # homogeneous using QH and I2.

    # (4) If any class is heterogeneous,
    # then the step 2 will be repeated with the rest of moderators



    # Then, step 2 will be repeated until we reach one of these three situations:
    #   (1) to reach a subset of effect sizes that
    # is homogeneous under the fixed effect model,
    # (2) to reach a subset of effect sizes that is still heterogeneous
    # but that it is too small to keep partitioning, or
    # (3) to reach a subset of effect sizes that is still heterogeneous
    # but the considered moderators are not able to
    # explain the heterogeneity in this subset.


    # (5) The partition process will stop when all classes are homogeneous,
    # when they are final due to small sample size,
    # or if none of the moderator explains the heterogeneity

  }



# Step 3: integrating effect sizes. ---------------------------------------


  list(Models, tree)


}
