####Main code###
#Load functions:
library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)


setwd("/Users/ratmir/MISE")


source("functions.R") # load all functions

#define cov before reading the functions
cov=0
sigma_sim=matrix(c(4, cov, cov,
               cov, 4, cov,
               cov, cov, 4), nrow=3, ncol=3)
#define model parameters before reading the functions
#beta_01 <<- 1; beta_02 <<- 1; beta_11 <<- 0.5; beta_12 <<- 0.5
#gamm <<- 5; con <<- 0

#Component-based:
#CASE1: Independence
source("integral_population.R")
source("integral_estimation.R")
source("shapley_int.R")
source("SE_vec_int.R")


#separate functions for DGP. (additive case)
g1 = function(X){ return( -sin(2*X[,1]) ) } #E(g1)=0. 
g2 = function(X){ return( cos(3*X[,2])  ) } #E(g2)=0 
g3 = function(X){ return( 0.5*X[,3] ) } #E(g3)=0. 
int = function(X){
  x1 = X[,1]
  x2 = X[,2]
  return( 2*cos(x1)*sin(2*x2)  ) 
}


#m_full_why = function(X){
#  x1=as.numeric(X[1])
#  x2=as.numeric(X[2])
#  x3=as.numeric(X[3])
#  return(
#     (-sin(2*x2) + cos(x3) )*(1/(1 + exp( -gamm*(x1 - con) ) ))
#         + (0.3*x2*x3 + cos(x2))*(1 - 1/(1 + exp( -gamm*(x1 - con) ) )) )
         
#}

#m_full_why = function(X){
#  x1=as.numeric(X[1])
#  x2=as.numeric(X[2])
#  x3=as.numeric(X[3])
#  return(beta_01*x2 + beta_02*x3 + beta_11*x2*(1/(1 + exp( -gamm*(x1 - con) ) ))
#         + beta_12*x3*(1/(1 + exp( -gamm*(x1 - con) ) )) )
#}

m_full_why = function(X){
  x1=as.numeric(X[1])
  x2=as.numeric(X[2])
  x3=as.numeric(X[3])
  return(-sin(2*x1) + cos(3*x2) + 0.5*x3 + 2*cos(x1)*sin(2*x2) )
}


true_model_list = list()
true_model_list[[1]] = m_x1
true_model_list[[2]] = m_x2
true_model_list[[3]] = m_x3 
true_model_list[[4]] = m_x1_x2
true_model_list[[5]] = m_x1_x3
true_model_list[[6]] = m_x2_x3 
true_model_list[[7]] = m_full_why 

#true_model_list = list()
#true_model_list[[1]] = g1
#true_model_list[[2]] = g2
#true_model_list[[3]] = g3 
#true_model_list[[4]] = function(X){return(g1(X) + g2(X)) }
#true_model_list[[5]] = function(X){return(g1(X) + g3(X))}
#true_model_list[[6]] = function(X){return(g2(X) + g3(X))}
#true_model_list[[7]] = function(X){return(g1(X) + g2(X) + g3(X)) }







l = -2; u = 2; N=1000; M=1
l_int = l; u_int = u
#ISE1 = rep(0, M) 
#ISE2 = rep(0, M) 
#ISE3 = rep(0, M) 
d = 3





#y=sapply(seq(-1,1, length.out=500), function(a)
#  {return(a ) # -sin(4a) macht sinn, cos(3a)
  
#})

#plot(y, x=seq(-1,1, length.out=500), type="l")

ISE_fct = function(m){
  
X<<-data.frame(mvrnorm(n=N, mu=c(0,0,0), Sigma=sigma_sim))

#DGP
Y <<- g1(X) + g2(X) + g3(X) + int(X) + rnorm(nrow(X), mean=0, sd=1)
#Y <<- g1(X) + g2(X) + g3(X) + int(X) + rt(n=nrow(X), df=10)



#All possible subsets
subs <<- subsets(X)

#Get model fits and sort them in a list
model_list <<- model_list_fct(subs=subs) # 75sek, 14 sek with tol = 0.1 and ftol = 0.1 and 1 multistart

while (sum(model_list[[7]]$bw[1:2]>10)>0
       | sum(model_list[[6]]$bw[1]>10)>0 
       | sum(model_list[[5]]$bw[1]>10)>0
       | sum(model_list[[4]]$bw[1:2]>10)>0
       | sum(model_list[[2]]$bw[1]>10)>0
       | sum(model_list[[1]]$bw[1]>10)>0)
  {
  X<<-data.frame(mvrnorm(n=N, mu=c(0,0,0), Sigma=sigma_sim))
  Y <<- g1(X) + g2(X) + g3(X) + rnorm(nrow(X))
  model_list <<- model_list_fct(subs=subs)
}



  


print(model_list[[7]])
# Component-based
ISE_res1=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1)
ISE_res2=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2)
ISE_res3=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3)

ISE1 = ISE_res1$integral 
ISE2 = ISE_res2$integral 
ISE3 = ISE_res3$integral 

#print(c(ISE1, ISE2, ISE3))


#system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1))[3]
#system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=1e-1, j=1))[3]

#system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2))[3]
#system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=1e-1, j=2))[3]

#system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3))[3]
#system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=1e-1, j=3))[3]




# Integral-based
ISE_res1_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1)# 45sek für 0.3, gleich 0.5, intern zusätzl 0.3 gibt 31
ISE_res2_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2)# 47sek für 0.3, gleich 0.5, intern zusätzl 0.3 gibt 34
ISE_res3_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3)# 100sek für 0.3, 15 sek 0.5, intern zusätzl 0.3 gibt 11

ISE1_int = ISE_res1_int$integral 
ISE2_int = ISE_res2_int$integral  
ISE3_int = ISE_res3_int$integral


#print(c(ISE1_int, ISE2_int, ISE3_int))

return(c(ISE1, ISE2, ISE1_int, ISE2_int))
#return(1)

}

collect=matrix(0, nrow=4, ncol=4)
for (m in 1:20){
  ISE_fct(m)
}

#N=200 (mit multi start, OHNE cfac.init=1, ftol=0.1)
#> colMeans(test)
#[1] 13.1529672  9.0987185  2.0647313 24.4285323  4.6165404  0.4832041
#> colMeans(collect)
#[1] 12.8018048  7.6503392  3.1848473 19.9955716  3.3656189  0.5903145

#N=500
#4.4098690 4.2162349 1.6463700 2.9109092 3.2249723 0.4562511
#N=1000
#2.5967004 2.5946939 0.8902474 2.0789921 2.2972707 0.4260488


library(parallel)
results_list = mclapply(1:8, ISE_fct, mc.cores=4)
system.time(mclapply(1:4, ISE_fct, mc.cores=4))[3]

results = matrix(unlist(results_list), byrow = FALSE, ncol=8)
system.time(sapply(1:10, ISE_fct))[3]




#comparison of computational integral specs + runtime check:
var=3
a=cubintegrate(f = SE_vec, lower = rep(l_int,d), upper = rep(u_int,d), method = "cuhre",
                      relTol = 1e-1, j=var) 
a.time = system.time(cubintegrate(f = SE_vec, lower = rep(l_int,d), upper = rep(u_int,d), method = "cuhre",
                                  relTol = 1e-1, j=var) )[3]
b=cubintegrate(f = SE_vec, lower = rep(l_int,d), upper = rep(u_int,d), method = "cuhre",
                      relTol = 1e-1, j=var, nVec = 128L) 
b.time = system.time(cubintegrate(f = SE_vec, lower = rep(l_int,d), upper = rep(u_int,d), method = "cuhre",
                                  relTol = 1e-1, j=var, nVec = 128L) )[3]
c=cubintegrate(f = SE_vec, lower = rep(l_int,d), upper = rep(u_int,d), method = "cuhre",
             relTol = 1e-2, j=var, nVec = 128L) 
c.time = system.time(cubintegrate(f = SE_vec, lower = rep(l_int,d), upper = rep(u_int,d), method = "cuhre",
                                  relTol = 1e-2, j=var, nVec = 128L) )[3]
de=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=1e-1, j=var)
de.time = system.time(hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=1e-1, j=var))[3]
e=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=1e-1, j=var, vectorInterface = TRUE)
e.time = system.time(hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=1e-1, j=var, vectorInterface = TRUE))[3]
# hcubature, hcubature, hcubature

a$integral
b$integral
c$integral
de$integral
e$integral

a.time
b.time
c.time
de.time
e.time














#Make a 3D plot to double check that shapley curve is correctly calculated.
# go in shapley_popul_vec(j, x_eval), where x_eval is 


x1_grid = seq(-2, 2, length.out=30) 
x2_grid = seq(-2, 2, length.out=30)


grid=t(expand.grid(x1_grid, x2_grid))
grid = rbind(grid , rep(0,ncol(grid)))
#evaluate on pop

shap_eval_est1 = matrix(0, nrow=ncol(grid), ncol=1)
shap_eval_est2 = matrix(0, nrow=ncol(grid), ncol=1)

for (i in 1:ncol(grid)){
  print(i)
  grid_col=as.numeric(grid[,i])
 shap_eval_est1[i] = shapley_int(j=1, grid_col) 
 shap_eval_est2[i] = shapley_int(j=2, grid_col) 
 
}

#shap_eval1=shapley_vec(j=1, grid) 
#shap_eval2=shapley_vec(j=2, grid)

#evaluate on est
shap_eval_est1=shapley_vec(j=1, grid) 
shap_eval_est2=shapley_vec(j=2, grid)

#shap_eval1=shapley_vec(j=1, grid) 
#shap_eval2=shapley_vec(j=2, grid)


shap_eval1 = matrix(0, nrow=ncol(grid), ncol=1)
shap_eval2 = matrix(0, nrow=ncol(grid), ncol=1)

for (i in 1:ncol(grid)){
  print(i)
  grid_col=as.numeric(grid[,i])
  shap_eval1[i] = shapley_popul(j=1, grid_col) 
  shap_eval2[i] = shapley_popul(j=2, grid_col) 
  
}

#shap_eval_est1=shapley_popul_vec(j=1, grid) 
#shap_eval_est2=shapley_popul_vec(j=2, grid)


shap_SE1 = SE_vec(j=1, grid)
surface1_SE=t(pracma::Reshape(shap_SE1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 

shap_SE2 = SE_vec(j=2, grid)
surface2_SE=t(pracma::Reshape(shap_SE2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 





surface1=t(pracma::Reshape(shap_eval1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface2=t(pracma::Reshape(shap_eval2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface1_est=t(pracma::Reshape(shap_eval_est1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface2_est=t(pracma::Reshape(shap_eval_est2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 



library(plotly)
par(mar = c(0, 0, 0, 2))
fig1 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1, type="surface") %>% hide_colorbar()
fig1 <- fig1 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap1')))

fig1_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1_est, type="surface") %>% hide_colorbar()
fig1_est <- fig1_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap1')))

fig2 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2, type="surface") %>% hide_colorbar()
fig2 <- fig2 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap2')))

fig2_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2_est, type="surface") %>% hide_colorbar()
fig2_est <- fig2 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap2')))

fig1_SE = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1_SE, type="surface") %>% hide_colorbar()
fig1_SE <- fig1_SE %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='SE shap1')))

fig2_SE = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2_SE, type="surface") %>% hide_colorbar()
fig2_SE <- fig2_SE %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='SE shap2')))



plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface1 , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,140)","rgb(0,3,140)"))) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface1_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))    ))

plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface2 , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,140)","rgb(0,3,140)"))) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface2_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))      ))


fig1_SE
fig2_SE




#setwd("/Users/ratmir/rdata")
#save(x1_grid, x2_grid, surface1, surface2, surface1_est, surface2_est, 
#     surface1_est_integ, surface2_est_integ, file = "data_GK.RData")








x=seq(from=-2,to=2,length.out=300)
y=-0.5^(sqrt(abs(x)))*cos(2*x) 
plot(y=(x^2 + x^3)/max(x^2 + x^3),x=x)
plot(y=x^3/max(x^3) , x=x, type="l")

plot(y= x, x=x)

#p=0, N(0,1), non-linear
#N=300
r1=c(9.276155, 12.148048,  1.662530, 12.948843, 14.600512,  0.441395)
r1_c= c(r1[1], r1[4], r1[2], r1[5], r1[3], r1[6])

#N=500
r2=c(5.824847, 6.740869, 0.881663, 7.573807, 7.649030, 0.247220)
r2_c= c(r2[1], r2[4], r2[2], r2[5], r2[3], r2[6])

#N=1000
r3=c(3.1361946, 4.0334151, 0.5123083, 4.0019116, 4.6224624, 0.1363247)
r3_c= c(r3[1], r3[4], r3[2], r3[5], r3[3], r3[6])

#N=2000
r4=c(1.95401555, 2.20016289, 0.31711195, 2.46146815, 2.57198453, 0.08280192)
r4_c= c(r4[1], r4[4], r4[2], r4[5], r4[3], r4[6])

tab1=rbind(r1_c, r2_c, r3_c, r4_c)
rownames(tab1)=c(300, 500, 1000, 2000)

#p=0.8, N(0,1), non-linear
#N=300
r1=c(11.186707, 11.699624,  3.205317, 12.409150, 14.375116,  1.882695)
r1_c= c(r1[1], r1[4], r1[2], r1[5], r1[3], r1[6])

#N=500
r2=c(8.004975,  8.104792,  2.275436,  7.874704, 10.102461,  1.309492)
r2_c= c(r2[1], r2[4], r2[2], r2[5], r2[3], r2[6])

#N=1000
r3=c(5.0807621, 4.6698189, 1.4877197, 4.7379844, 5.9228996, 0.7520788)
r3_c= c(r3[1], r3[4], r3[2], r3[5], r3[3], r3[6])

#N=2000
r4=c(3.3634134, 2.9993775, 1.0178001, 2.8826786, 3.6834617, 0.4997978)
r4_c= c(r4[1], r4[4], r4[2], r4[5], r4[3], r4[6])

tab2=rbind(r1_c, r2_c, r3_c, r4_c)
rownames(tab2)=c(300, 500, 1000, 2000)

#p=0, N(0,1), additive
#N=300
r1=c(8.463692,  8.683261,  3.227080, 13.020629,  9.457465,  0.375672)
r1_c= c(r1[1], r1[4], r1[2], r1[5], r1[3], r1[6])

#N=500
r2=c(5.1849977, 5.8812080, 2.0025507, 5.9693650, 6.5782766, 0.2080727)
r2_c= c(r2[1], r2[4], r2[2], r2[5], r2[3], r2[6])

#N=1000
r3=c(3.0860187, 3.4456492, 1.1003710, 3.2108246, 3.8596743, 0.1173032)
r3_c= c(r3[1], r3[4], r3[2], r3[5], r3[3], r3[6])

#N=2000
r4=c(1.829050, 2.117116, 0.66449743, 1.973804, 2.338630, 0.06675618)
r4_c= c(r4[1], r4[4], r4[2], r4[5], r4[3], r4[6])
tab3=rbind(r1_c, r2_c, r3_c, r4_c)
rownames(tab3)=c(300, 500, 1000, 2000)

#p=0.8, N(0,1), additive
#N=300
r1=c(7.223271, 8.702114, 2.530023, 7.926072, 9.377656, 1.289431)
r1_c= c(r1[1], r1[4], r1[2], r1[5], r1[3], r1[6])

#N=500
r2=c(5.1905900, 5.8234125, 1.8474942, 5.3735086, 6.5978481, 0.8754296)
r2_c= c(r2[1], r2[4], r2[2], r2[5], r2[3], r2[6])

#N=1000
r3=c(3.3579009, 3.5604188, 1.4330502, 3.3873674, 3.8363271, 0.5579934)
r3_c= c(r3[1], r3[4], r3[2], r3[5], r3[3], r3[6])

#N=2000
r4=c(2.2228183, 2.1698058, 0.8782387, 2.1537379, 2.3270589, 0.3188351)
r4_c= c(r4[1], r4[4], r4[2], r4[5], r4[3], r4[6])
tab4=rbind(r1_c, r2_c, r3_c, r4_c)
rownames(tab4)=c(300, 500, 1000, 2000)

##################################################
#p=0, t(10), non-linear
#N=300
r1=c(10.4122919, 13.6588072,  2.3020329, 14.9259994, 17.9368541,  0.4759117)
r1_c= c(r1[1], r1[4], r1[2], r1[5], r1[3], r1[6])

#N=500
r2=c(6.4073562, 7.4761741, 0.9781091, 8.8447628, 8.8501185, 0.2741349)
r2_c= c(r2[1], r2[4], r2[2], r2[5], r2[3], r2[6])

#N=1000
r3=c(3.4966011, 4.3742736, 0.5541970, 4.6988047, 5.1534854, 0.1483494)
r3_c= c(r3[1], r3[4], r3[2], r3[5], r3[3], r3[6])

#N=2000
r4=c(2.17869482, 2.42245387, 0.31046600, 2.82573790, 3.01564699, 0.09125454)
r4_c= c(r4[1], r4[4], r4[2], r4[5], r4[3], r4[6])

tab5=rbind(r1_c, r2_c, r3_c, r4_c)
rownames(tab5)=c(300, 500, 1000, 2000)

#p=0.8, t(10), non-linear
#N=300
r1=c(10.437936, 12.104519,  3.402987, 12.032228, 14.822777,  2.134425)
r1_c= c(r1[1], r1[4], r1[2], r1[5], r1[3], r1[6])

#N=500
r2=c(8.475416,  8.399251,  2.150528,  9.269826, 10.392018, 1.453637)
r2_c= c(r2[1], r2[4], r2[2], r2[5], r2[3], r2[6])

#N=1000
r3=c(5.4705635, 5.5356351, 1.6266648, 5.3557430, 6.4949406, 0.8939049)
r3_c= c(r3[1], r3[4], r3[2], r3[5], r3[3], r3[6])

#N=2000
r4=c(3.4552666, 3.2961019, 1.0701620, 3.2165842, 3.9595917, 0.5411341)
r4_c= c(r4[1], r4[4], r4[2], r4[5], r4[3], r4[6])

tab6=rbind(r1_c, r2_c, r3_c, r4_c)
rownames(tab6)=c(300, 500, 1000, 2000)

#p=0, t(10), additive
#N=300
r1=c(10.3402425,  8.3874048,  2.8487331, 19.1584805,  6.6665210,  0.3298615)
r1_c= c(r1[1], r1[4], r1[2], r1[5], r1[3], r1[6])

#N=500
r2=c(5.5083132, 6.0815622, 1.9104747, 6.4888730, 7.6198076, 0.2304462)
r2_c= c(r2[1], r2[4], r2[2], r2[5], r2[3], r2[6])

#N=1000
r3=c(3.3808654, 3.7635950, 1.1329254, 3.8577505, 4.5410704, 0.1248227)
r3_c= c(r3[1], r3[4], r3[2], r3[5], r3[3], r3[6])

#N=2000
r4=c(1.94502188, 2.35287129, 0.64988993, 2.23288318, 2.74880918, 0.07828728)
r4_c= c(r4[1], r4[4], r4[2], r4[5], r4[3], r4[6])
tab7=rbind(r1_c, r2_c, r3_c, r4_c)
rownames(tab7)=c(300, 500, 1000, 2000)

#p=0.8, t(10), additive
#N=300
r1=c(7.675637,  8.890430,  3.109195,  8.559716, 10.191730,  1.520296)
r1_c= c(r1[1], r1[4], r1[2], r1[5], r1[3], r1[6])

#N=500
r2=c(5.459700, 5.981454, 2.067453, 5.778860, 7.118244, 1.018030)
r2_c= c(r2[1], r2[4], r2[2], r2[5], r2[3], r2[6])

#N=1000
r3=c(3.5469635, 3.6700762, 1.3750320, 3.7128736, 4.1925309, 0.6367359)
r3_c= c(r3[1], r3[4], r3[2], r3[5], r3[3], r3[6])

#N=2000
r4=c(2.5793250, 2.5173634, 0.8500102, 2.4838458, 2.7687481, 0.3698739)
r4_c= c(r4[1], r4[4], r4[2], r4[5], r4[3], r4[6])
tab8=rbind(r1_c, r2_c, r3_c, r4_c)
rownames(tab8)=c(300, 500, 1000, 2000)







