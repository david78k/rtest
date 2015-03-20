library(recommenderlab)
library(ggplot2)

data(MovieLense)
MovieLense

#png(filename = "rawratings.png")
image(sample(MovieLense, 500), main = "Raw ratings")

#png(filename = "histogram.png")
qplot(getRatings(MovieLense), binwidth = 1, main = "Histogram of ratings", xlab = "Rating")
summary(getRatings(MovieLense)) # skewed to the right

qplot(rowCounts(MovieLense), bindwidth = 10, 
	main = "Movies Rated on average",
	xlab = "# of users",
	ylab = "# of movies rated")
