library("markovchainTaeSeungKang")
#library("markovchain")

sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
#sequence <- data.frame(t(sequence))
#microbenchmark(
#  markovchainFit(data = sequence)
  #markovchainFit(data = sequence, method="laplace", laplacian=0.1),
  #markovchainFit(data = sequence, method="bootstrap"),
  #mcfit(data = sequence, method="bootstrap"),
  #markovchainFit(data = sequence, byrow=FALSE)#,

  markovchainFit_cpp(sequence)
  #markovchainFit_cpp(sequence, "laplace", laplacian=0.1)
  #markovchainFit_cpp(sequence, "bootstrap")
  #markovchainFit_cpp(sequence, byrow=FALSE)
#)

