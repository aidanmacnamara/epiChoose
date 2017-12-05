# example showing the difference between mean and correlation

x = data.frame(
  s1=rnorm(20,0,1),
  s2=rnorm(20,0,1)
)

# x = data.frame(
#   s1=rnorm(20,0,1),
#   s2=rnorm(20,2,1)
# )

# s1 = rnorm(20)
# s2 = s1+2+rnorm(20,0,0.5)
# x = data.frame(s1, s2)

# s1 = rnorm(20)
# s2 = s1+rnorm(20,0,0.5)
# x = data.frame(s1, s2)

p1 = ggplot(x, aes(s1, s2)) + geom_point() + theme_thesis() + xlab("") + ylab("")
p2 = ggplot(melt(x), aes(variable, value)) + theme_thesis() + geom_boxplot() + xlab("") + ylab("")
plot_grid(p1, p2)

t.test(x$s1, x$s2, paired=T)$p.value
cor.test(x$s1, x$s2)$p.value
