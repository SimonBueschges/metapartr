## This is the global script, which loads all necessary sources,
## sets the options for the application and loads all the data
## which can be loaded before the start of any session.
## It also prepares this loaded data for later usage.
##
##

##### Clear Workspace ##########################################################

rm(list = ls())

##### options ##################################################################

##### sources ##################################################################

##### > libraries ##############################################################

# Loading the libraries

##### >> Data treatment ########################################################
library(data.table)
library(magrittr)
library(stringr)
library(openxlsx)

##### >> Assertions ############################################################
library(checkmate)


##### >> Visualisation #########################################################

library(ggplot2)
library(ggcorrplot)
library(DiagrammeR)

##### >> Report Creation #######################################################
library(rmarkdown)
library(rticles)
library(xtable)
library(kableExtra)


##### >> Statistics ############################################################

library(psych)
library(metafor)
library(meta)


library(metapartr)
##### > functions ##############################################################

source("R/meta_partition.r", encoding = "UTF-8")

source("R/heterogeneity_candidates.r", encoding = "UTF-8")

##### data import ##############################################################

# Assure Existince
# assert_file_exists("./Data/Extracted Data.xlsx")

# Data and Tables
Studies <-
  read.xlsx("./Data/DataForMethod.xlsx") %>%
  as.data.table() %>%
  setorder(
    "Study",
    "ESr",
    "n",
    "Country",
    "Species",
    "Taxa",
    "Study.Type",
    "Sampling.Effort",
    "Patch.Area"
  )


# I see 229 Species, but they say 220?
Species.Info <-
  read.xlsx("./Data/MetaInfoForMethodData.xlsx") %>%
  as.data.table() %>%
  .[, `Home.Range.(ha)` := as.numeric(`Home.Range.(ha)`)] %>%
  .[, `Length.(cm)` := as.numeric(`Length.(cm)`)] %>%
  .[, `Repro` := as.numeric(`Repro`)]


Method.Study.Data <-
  merge.data.table(Studies,
    Species.Info,
    by = c("Study", "Taxa", "Species")
  ) %>%
  setorder(
    "Study",
    "Taxa",
    "Species",
    "Order",
    "Country",
    "ESr",
    "n",
    "Study.Type",
    "Sampling.Effort",
    "Patch.Area",
    "Repro",
    "Home.Range.(ha)",
    "Length.(cm)"
  ) %>%
  na.omit()

##### data preparation #########################################################
#
# Method.Study.Data %>%
#   write.xlsx("data/Quesnelle.Species.Sensitivity.xlsx")



##### exploration ##############################################################

##### modeling #################################################################
#
Model <-
  metacor(ESr, n, data = Method.Study.Data, studlab = Study)

formula <-
  TE ~
  Study.Type +
  Sampling.Effort +
  Patch.Area +
  `Home.Range.(ha)` +
  `Length.(cm)` +
  Repro +
  Taxa


Metapart <-
  meta_partition(
    model = Model,
    moderators = c(
      "Study.Type", "Sampling.Effort",
      "Patch.Area",
      "Home.Range.(ha)", "Length.(cm)", "Repro",
      "Taxa"
    ),
    n.to.small = 8,
    auto = TRUE,
    c = 1,
    control = rpart.control(
      xval = 10,
      minbucket = 3,
      minsplit = 6,
      cp = 0.0001
    )
  )

saveRDS(Metapart, "data/Metapart.rds")

# library(Rcpp)
# FEtree <- FEmrt(formula, vi = seTE, data = cbind(Method.Study.Data,"TE"= Model$TE,"seTE" = Model$seTE), c = 0)
# print(FEtree)
# summary(FEtree)
# plot(FEtree)
#
# svg(filename = "graphics/Method.Study.Forest.svg",width =  12, height = 96)
# Model %>%
#   forest.meta()
# dev.off()
#
# Method.Study.Data %>%
#   names %>%
#   .[!(. %in% c(
#     "Study",
#     "ESr",
#     "n",
#     "Repro",
#     "Mass",
#     "Length",
#     "Home.Range"
#   ))] %>%
#   purrr::walk(
#     function(x){
#
#       Model <-
#         metacor(ESr,
#                 n ,
#                 data = Method.Study.Data,
#                 studlab = Study,
#                 byvar = get(x)
#                 )
#
#
#       svg(filename = paste0("graphics/Method.Study.Forest.", x, ".svg"),
#           width =  12,
#           height = 96)
#       Model %>%
#         forest.meta()
#       dev.off()
#     }
#   )


# a <- Metapart %>%
#    sapply(function(x){
#      if(is.null(x$Terminal.Node)){FALSE
#        }else{
#          x$Terminal.Node}
#      })
# A <- Metapart[a] %>% lapply(function(x){x[["Model"]]$subset})
# b <- vector()
# for(i in 1:27){
#   b[A[[i]]] <- i
# }
#
# svg(filename = paste0("graphics/Method.Study.Forest.method.svg"),
#     width =  12,
#     height = 96)
# update.meta(Model, byvar = a) %>%
#   forest.meta(sortvar = TE)
# dev.off()
