
#################################################
# BELOW ARE FUNCTION VERSIONS FOR SAMPLE NICHES #
#################################################

individual_contribution <- function(x, w = NULL) {
  # x is the individual parameter data, columns are mu and sigma (note that it is the sd not the variance),
  # but the returned population sigma2 is the variance
  n = nrow(x)
  if(!is.null(w)) n = w
  mu_i = as.numeric(x[,"mu"])
  sigma_i = as.numeric(x[,"sigma"])
  if("skew" %in% colnames(x)){
    skew_i = as.numeric(x[,"skew"])
  }
  mu_i = mu_i - mean(mu_i)
  mu = mean(mu_i)
  marginality_sigma2 = (mu_i^2 - mu^2)/n
  specialization_sigma2 = (sigma_i^2)/n
  sigma2 = sum(marginality_sigma2) + sum(specialization_sigma2)
  marginality_skew = (mu_i^3 - mu^3)/(n*sigma2^(3/2))
  specialization_skew = (3*mu_i*(sigma_i^2)-3*mu*sigma2)/(n*sigma2^(3/2))
  if("skew" %in% colnames(x)){
    indiv_skew = (skew_i*(sigma_i^3))/(n*sigma2^(3/2))
    specialization_skew = specialization_skew + indiv_skew
  }
  skew = sum(marginality_skew) + sum(specialization_skew)
  
  output = cbind(mu_pop = mu,marginality_sigma2 = marginality_sigma2,
                 specialization_sigma2 = specialization_sigma2,
                 sigma2_pop = sigma2,
                 marginality_skew = marginality_skew,
                 specialization_skew = specialization_skew,
                 skew_pop = skew)
  if("skew" %in% colnames(x)) output = cbind(output,indiv_skew)
  return(output)
}

# computes the individual densities of all individuals along w/ the mixture dist'n of all indivs
mixture_densities0 <- function(individual_niche_params,
                               x = seq(-10, 10, 0.1)) {
  indiv_dens = apply(individual_niche_params, 1, function(y){
    dnorm(x,y[1],y[2])
  })
  mix_den = apply(indiv_dens, 1, mean)
}

indiv_densities0 <- function(individual_niche_params,
                             x = seq(-10, 10, 0.1)) {
  indiv_dens = apply(individual_niche_params, 1, function(y){
    dnorm(x,y[1],y[2])
  })
}


#################################################
# Conceptual figures #
# nested
#pdf("Fig1_conceptual_examples_2022_03_29.pdf",width=8,height=8)
pdf("Fig1_conceptual_examples_2023_06_25.pdf",width=6,height=6)
layout.matrix <- matrix(c(1:9), nrow = 3, ncol = 3, byrow = F)
layout(mat = layout.matrix) # Widths of the two columns
par(mar=c(3.5,2,0.5,2))
x = seq(-10, 10, 0.1)+24
# indivs = cbind(mu = c(-2,-1,0,1,2,3)+24,
#                sigma = c(4,5,1,3,6,2)^0.5)
# mix = mixture_densities0(indivs,x)
# indiv_dens = indiv_densities0(indivs,x)
# indiv.cont = individual_contribution(indivs)
# colors = c("dodger blue","light green","orange","cyan",'blue violet',"red")
# for (i in 1:ncol(indiv_dens)) {
#   if(i==1) plot(indiv_dens[,i]~x,type='l',lwd=2,ylim=range(0,0.5),col=colors[i],main="",xlab='',ylab='')
#   lines(indiv_dens[,i]~x,lwd=2,col=colors[i])
# }
# lines(mix~x,lwd=3,col=1)
# par(mar=c(2,2,2,2))
# #legend("topleft",c(paste0("INDV_",1:6),"POP"),lwd=1.5,col=c(2:7,1),cex=0.8)
# plot(indiv.cont[,c("marginality_sigma2","specialization_sigma2")],xlab="marginality_variance",ylab="specialization_variance",
#      col=colors,cex=2.5,pch=16,main='')
# plot(indiv.cont[,c("marginality_skew","specialization_skew")],xlab="marginality_skewness",ylab="specialization_skewness",
#      col=colors,cex=2.5,pch=16,main='')
# abline(h=0,v=0,lty="dashed",col="grey")

par(mar=c(3.5,2,0.5,2))
indivs = cbind(mu = (-2):3+24,sigma = 2)
mix = mixture_densities0(indivs,x)
indiv_dens = indiv_densities0(indivs,x)
indiv.cont = individual_contribution(indivs)
for (i in 1:ncol(indiv_dens)) {
  if(i==1) plot(indiv_dens[,i]~x,type='l',lwd=2,ylim=range(0,0.3),ylab="",col=colors[i],main=" ",xlab='',bty = 'l')
  lines(indiv_dens[,i]~x,lwd=2,col=colors[i])
}
lines(mix~x,lwd=3,col=1)
#legend("topleft",c(paste0("INDV_",1:6),"POP"),lwd=1.5,col=c(2:7,1),cex=0.8)
par(mar=c(2,2,2,2))
indiv.cont[,c("marginality_sigma2")] = indiv.cont[,c("marginality_sigma2")] + c(rep(0.025,3),rep(0,3))
plot(indiv.cont[,c("marginality_sigma2","specialization_sigma2")],xlab="marginality_variance",ylab="specialization_variance",
     col=colors,cex=2.5,pch=16,main='',bty='l')
plot(indiv.cont[,c("marginality_skew","specialization_skew")],xlab="marginality_skewness",ylab="specialization_skewness",
     col=colors,cex=2.5,pch=16,main='',bty='l')
abline(h=0,v=0,lty="dashed",col="grey")

par(mar=c(3.5,2,0.5,2))
indivs = cbind(mu = rep(c(-1,2),each=3)+24,
               sigma = rep(c(3,2,1),2))
mix = mixture_densities0(indivs,x)
indiv_dens = indiv_densities0(indivs,x)
indiv.cont = individual_contribution(indivs)
for (i in 1:ncol(indiv_dens)) {
  if(i==1) plot(indiv_dens[,i]~x,type='l',lwd=2,ylim=range(0,0.5),ylab="",col=colors[i],main="",xlab='',bty='l')
  lines(indiv_dens[,i]~x,lwd=2,col=colors[i])
}
lines(mix~x,lwd=3,col=1)
#legend("topleft",paste0("INDV_",1:6),lwd=1.5,col=2:7)
par(mar=c(2,2,2,2))
indiv.cont[,c("specialization_sigma2")] = indiv.cont[,c("specialization_sigma2")] + c(rep(0.05,3),rep(0,3))
plot(indiv.cont[,c("marginality_sigma2","specialization_sigma2")],xlab="marginality_variance",ylab="specialization_variance",
     col=colors,cex=2.5,pch=16,main='',bty='l')
plot(indiv.cont[,c("marginality_skew","specialization_skew")],xlab="marginality_skewness",ylab="specialization_skewness",
     col=colors,cex=2.5,pch=16,main='',bty='l')
abline(h=0,v=0,lty="dashed",col="grey")

par(mar=c(3.5,2,0.5,2))
indivs = cbind(mu = (-2):3+24,
               sigma = (6:1)^0.5)
mix = mixture_densities0(indivs,x)
indiv_dens = indiv_densities0(indivs,x)
indiv.cont = individual_contribution(indivs)
for (i in 1:ncol(indiv_dens)) {
  if(i==1) plot(indiv_dens[,i]~x,type='l',lwd=2,ylim=range(0,0.5),ylab="",col=colors[i],main="",xlab='',bty='l')
  lines(indiv_dens[,i]~x,lwd=2,col=colors[i])
}
lines(mix~x,lwd=3,col=1)
#legend("topleft",paste0("INDV_",1:6),lwd=1.5,col=2:7)
par(mar=c(2,2,2,2))
plot(indiv.cont[,c("marginality_sigma2","specialization_sigma2")],xlab="marginality_variance",ylab="specialization_variance",
     col=colors,cex=2.5,pch=16,main='',bty='l')
plot(indiv.cont[,c("marginality_skew","specialization_skew")],xlab="marginality_variance",ylab="specialization_variance",
     col=colors,cex=2.5,pch=16,main='',bty='l')
abline(h=0,v=0,lty="dashed",col="grey")
dev.off()
