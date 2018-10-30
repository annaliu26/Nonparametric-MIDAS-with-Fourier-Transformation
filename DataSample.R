load('I:/My Drive/Research/Nonparametric estimation/Simulation Results/DGP-sample.RData')
data = DATA.e[[1]]

R = 15
y1 = as.matrix(data$DATA1[[1]][,1]); X1 = as.matrix(data$DATA1[[1]][,-1])
y4 = as.matrix(data$DATA4[[1]][,1]); X4 = as.matrix(data$DATA4[[1]][,-1])

for (r in 2:R)
{
  y11 = as.matrix(data$DATA1[[r]][,1]); X11 = as.matrix(data$DATA1[[r]][,-1])
  y44 = as.matrix(data$DATA4[[r]][,1]); X44 = as.matrix(data$DATA4[[r]][,-1])
  y1 = rbind(y1,y11); X1 = rbind(X1,X11)
  y4 = rbind(y4,y44); X4 = rbind(X4,X44)
}
y = rbind(y1,y4); X = rbind(X1,X4)
table = cbind(y,X)

write.csv(table, 'I:/My Drive/Research/Nonparametric estimation/data.csv',row.names=FALSE)


