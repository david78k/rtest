library("markovchain")
library("Rcpp")

showClass("markovchain")
transmat <- matrix(c(0.4, 0.6,.3, .7), nrow = 2, byrow = TRUE)
simplemc <- new("markovchain", states=c("a", "b"),
		transitionMatrix = transmat,
		name = "simplemc")

simplemc^4

steadyStates(simplemc)
absorbingStates(simplemc)
simplemc[2, 1]
t(simplemc)
is.irreducible(simplemc)
conditionalDistribution(simplemc, "b")

sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcfit <- markovchainFit(data = sequence)

predict(mcfit$estimate, newdata = "b", n.ahead = 3)
mymc <- as(transmat, "markovchain")

summary(simplemc)

#plot(simplemc)

