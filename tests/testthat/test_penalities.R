context("Test Penalties")

# set the seed and generate a tree 
set.seed(42)
neg=100/365
off.r=5
w.shape=10
w.scale=0.1
pi=0.25
dateT=2008
simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
                         w.scale=w.scale,dateStartOutbreak=2005,dateT=dateT)

case.info <- simulateExposure(simu, propCases = 1)
contact.info <- simulateContact(simu, propCases = 1)

test_that("test correct simulated exposure", {
  epiData <- list(case.info = case.info)
  
  # should be no pen 
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),0)
  
  # should be no pen 
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = F)$penalties),0)
})

test_that("test correct simulated contact", {
  epiData <- list(contact.info = contact.info)
  
  # should be no pen 
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),0)
  
  # should be no pen 
  expect_equal(sum(unlist(epiPenTTree(extractTTree(simu), epiData, penaltyInfo = F))),0)
  
})

test_that("test correct exposure penalty", {
  epiData <- list(case.info = case.info)
  
  # offset the correct times enough so that they are all broken (+10 years should do)
  epiData$case.info[,1:2] <- epiData$case.info[,1:2] + 10
  
  # should be all pen  
  expect_equal(sum(epiPenTTree(extractTTree(simu), epiData = epiData, penaltyInfo = T)$penalties),nrow(case.info))
  
  # should be all pen 
  expect_equal(sum(unlist(epiPenTTree(extractTTree(simu), epiData = epiData, penaltyInfo = F))),nrow(case.info))
})

test_that("test correct contact penalty", {
  epiData <- list(contact.info = contact.info)
  
  # offset the correct times enough so that they are all broken (+10 years should do)
  epiData$contact.info[,3:4] <- epiData$contact.info[,3:4] + 10
  
  # should be all pen  
  expect_equal(sum(epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),nrow(contact.info))
  
  # should be all pen 
  expect_equal(sum(unlist(epiPenTTree(extractTTree(simu), epiData, penaltyInfo = F))),nrow(contact.info))
}) 


test_that("test correct location penalty", {
  # give every case a unique locaiton so that every present direct transmission is penalized 
  
  
}) 


# write tests to check that correct values are being placed into the info about who is breaking pen 



# test that just tracking does not change result 
# test_that("test track but no penalty", {
#   neg=100/365
#   off.r=5
#   w.shape=10
#   w.scale=0.1
#   pi=0.25
#   dateT=2008
#   set.seed(0)
#   record<-inferTTree(ptreeFromPhylo(rtree(5), dateLastSample = 2007.7),mcmcIterations=100,w.shape=w.shape,w.scale=w.scale,dateT=dateT)
#   
#   epiData <- list(case.info = simulateExposure(record[[100]]$ctree))
#   set.seed(0)
#   record2 <- inferTTree(ptreeFromPhylo(rtree(5), dateLastSample = 2007.7),mcmcIterations=100,w.shape=w.shape,w.scale=w.scale,dateT=dateT,epiData = epiData, penalize = F)
#   
#   expect_identical(record,record2)
#   
#   set.seed(0)
#   record2 <- inferTTree(ptreeFromPhylo(rtree(5), dateLastSample = 2007.7),mcmcIterations=100,w.shape=w.shape,w.scale=w.scale,dateT=dateT,epiData = epiData, penalize = F, trackPenalty = T)
#   for(i in 1:length(record2)){
#     record2[[i]]$penalty.information <- NULL
#     record2[[i]]$penalty.exposure <- NULL
#     record2[[i]]$penalty.contact <- NULL
#     record2[[i]]$penalty.location <- NULL
#   }
#   
#   expect_identical(record,record2)
# })

# test getting specific events wrong with epi data and 

