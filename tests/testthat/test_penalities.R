context("Test Penalties")

# set the seed and generate a tree to use in some tests 
set.seed(42)
neg=100/365
off.r=5
w.shape=10
w.scale=0.1
pi=0.25
dateT=2008
simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
                         w.scale=w.scale,dateStartOutbreak=2005,dateT=dateT)

exposure <- simulateExposure(simu, propCases = 1)
contact <- simulateContact(simu, propCases = 1)

test_that("test correct simulated exposure", {
  epiData <- list(exposure = exposure)
  
  # should be no pen 
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),0)
  
  # should be no pen 
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = F)),0)
})

test_that("test correct simulated contact", {
  epiData <- list(contact = contact)
  
  # should be no pen 
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),0)
  
  # should be no pen 
  expect_equal(sum(unlist(epiPenTTree(extractTTree(simu), epiData, penaltyInfo = F))),0)
  
})

test_that("test correct exposure penalty", {
  epiData <- list(exposure = exposure)
  
  # offset the correct times enough so that they are all broken (+10 years should do)
  epiData$exposure[,1:2] <- epiData$exposure[,1:2] + 10
  
  # should be all pen  
  expect_equal(sum(epiPenTTree(extractTTree(simu), epiData = epiData, penaltyInfo = T)$penalties),nrow(exposure))
  
  # should be all pen 
  expect_equal(sum(unlist(epiPenTTree(extractTTree(simu), epiData = epiData, penaltyInfo = F))),nrow(exposure))
})

test_that("test correct contact penalty", {
  epiData <- list(contact = contact)
  
  # offset the correct times enough so that they are all broken (+10 years should do)
  epiData$contact[,3:4] <- epiData$contact[,3:4] + 10
  
  # should be all pen  
  expect_equal(sum(epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),nrow(contact))
  
  # should be all pen 
  expect_equal(sum(unlist(epiPenTTree(extractTTree(simu), epiData, penaltyInfo = F))),nrow(contact))
}) 


test_that("test correct location penalty", {
  # give every case a unique locaiton so that every present direct transmission is penalized 
  
  nsam <- length(simu$nam)
  
  location <- as.character(1:nsam)
  names(location) <- simu$nam
  
  epiData <- list(location = location)
  
  ttree <- extractTTree(simu)
  
  wiw <- cbind(ttree$ttree[,3],1:length(ttree$ttree[,1]))
  
  wiw <- wiw[wiw[,2] != 0 & wiw[,1] != 0,]
  
  # filter to be only between sampled cases 
  sampledCaseRows <- wiw[,2] <= nsam & wiw[,1] <= nsam 
  
  # everyone has unique location shoudl all sam to sam should be pen 
  expect_equal(sum(unlist(epiPenTTree(extractTTree(simu), epiData, penaltyInfo = F))),sum(sampledCaseRows))
  
}) 

test_that("test track but no penalty", {
  
  ptree <- extractPTree(simu)
  
  neg=100/365
  off.r=5
  w.shape=10
  w.scale=0.1
  pi=0.25
  dateT=2008
  set.seed(0)
  record<-inferTTree(ptree,mcmcIterations=50,w.shape=w.shape,w.scale=w.scale,dateT=dateT)

  epiData <- list(exposure = simulateExposure(record[[50]]$ctree))
  set.seed(0)
  record_no_pen <- inferTTree(ptree,mcmcIterations=50,w.shape=w.shape,w.scale=w.scale,dateT=dateT,epiData = epiData, penalize = F)

  # test that just tracking does not change result
  expect_identical(record,record_no_pen)

  set.seed(0)
  record_track <- inferTTree(ptree,mcmcIterations=50,w.shape=w.shape,w.scale=w.scale,dateT=dateT,epiData = epiData, penalize = F, trackPenalty = T)
  
  # throw away the tracked info should be the same otherwise 
  for(i in 1:length(record_track)){
    record_track[[i]]$penalty.info <- NULL
    record_track[[i]]$penalty.exposure <- NULL
    record_track[[i]]$penalty.contact <- NULL
    record_track[[i]]$penalty.location <- NULL
  }

  expect_identical(record,record_track)
})

test_that("manual checks determined on fixed tree", {
  set.seed(3324234)
  simu <- simulateOutbreak(pi = 1)
  
  exposure <- data.frame(row.names = c(1,2,5), start = c(0,2000,2010), end = c(10,2005,2015))
  
  contact <- data.frame(case.A = c(2,5,1), case.B = c(4,1,3), start = c(2010,NA,2000), end = c(2011,NA,2005))
  
  location <- list("1" = c("a","b"),"3" = c("c","d"), "2"= "d","4" = "d")
  
  epi <- list(exposure = exposure, contact = contact, location = location)
  
  res <- epiPenTTree(extractTTree(simu), epiData = epi, penaltyInfo = T)
  
  # check quant. 
  expect_equal(unname(res$penalties),c(3,2,2))
  
  # exposure 
  expect_equal(unlist(unname(res$penalty.info$exposurePen[1,])),c("A",1,2,2001.14760960215,10,2000))
  
  # contact 
  expect_equal(round(unlist(unname(res$penalty.info$contactPen[1,]))),round(c(TRUE,2,4,2010,2011,2003.844)))
  
  # location 
  expect_equal(unlist(unname(res$penalty.info$locationPen[1,3])),c("a","b"))
  expect_equal(unlist(unname(res$penalty.info$locationPen[1,4])),c("d"))
  expect_equal(as.vector(unname(unlist(res$penalty.info$locationPen[1,1:2]))),c("1","2"))

})

test_that("test with NA values", {
  
  location <- as.character(1:(length(simu$nam)-5))
  names(location) <- simu$nam[1:(length(simu$nam)-5)]
  
  epiData <- list(exposure = exposure, contact = contact, location = location)
  
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),3)
  
  epiData$exposure[1,] <- NA
  epiData$exposure[2,1] <- NA
  epiData$exposure[3,2] <- NA
  
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),3)
  
  epiData$contact[1,3:4] <- NA
  epiData$contact[2,4] <- NA
  
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),3)
  
  epiData$exposure <- epiData$exposure + 100
  expect_equal(sum( epiPenTTree(extractTTree(simu), epiData, penaltyInfo = T)$penalties),21)
  
}) 
