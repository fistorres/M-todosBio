library("deSolve")

# Modelo CF sem tratamentos
CF_meso<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dB <- ((B^0.8) * (M^1.2) - (0.8*(B^1.2)) * (P^0.1))*M^(-2.4)
    dP <- 20*(B^0.2) * (D^0.4) * (A^(-0.1))-(50*P^0.5)
    dA <- 15*(B^0.1) * (P^0.1)-(12*A^0.5)
    dH <- 500-40*(H^0.3) *  (P^0.2) * (A^(- 0.2))
    dD <- (40*(H^0.3) * (P^0.2) * (A^( -0.2)))-5*D 
    dM <- 0.16*(2^CF) * (B^0.25) - 0.15*M^(2-CF)
    # return the rate of change
    list(c(dB,dP,dA,dH,dD,dM))
  }) # end with(as.list ...
}

# perturbações no valor de bactérias (infeções)
reset_Bac <- data.frame( 
  var = 'B', 
  time = c(20,160,300) , 
  value = c(3000,5000,10000),
  method = 'add')

# Quando CF=0, paciente saudável
parameters0 <- c(CF=0)
# valores iniciais de B,P,A,H,D,M (steady state) quando CF = 0
state0 <-  c(B=1.383617,P=6.091615,A=2.393149,H=2431.41,D=100,M=1.075577)
times <- seq(0, 400, by = 1)

out0 <- ode(y= state0,times = times, func = CF_meso,parms= parameters0, events = list(data=reset_Bac),method='impAdams')
#plot(out0)


Di <- 100
Hi = 2431.41
#lung function
LF <- 100*(out0[,'H']/out0[,'D'])/(Hi/Di)
#Percentagem de células H (healthy)
PH <- out0[,'H']*100/Hi

# recriação da figura 4.
par(mfrow=c(2,2))

plot(x=out0[,'time'],y=out0[,'D'],col='black',xlab='time',ylab="",type='l',ylim = c(65,120) )
lines(LF, col='red',lwd=1)
lines(PH, col='blue',lwd=1)
legend("bottomleft", legend=c('D','LF','PH'),
       col=c('black','red','blue'),horiz= TRUE,lty=1,cex = 0.74)

plot(x=out0[,'time'],y=out0[,'B'],col='black',xlab='time',ylab = 'B',type='l')

plot(x=out0[,'time'],y=out0[,'A'],col='black',xlab='time',ylab="",type='l',ylim=c(0,16))
lines(out0[,'M'], col='red',lwd=1)
legend('topleft', legend=c('A','M'),
       col=c('black','red'),lty=1,horiz=TRUE,cex=0.8)

plot(x=out0[,'time'],y=out0[,'P'],col='black',xlab='time',ylab = 'P',type='l')


# Quando CF=1, paciente doente.
parameters <- c(CF=1)
# valores iniciais de B,P,A,H,D,M (steady state) quando CF = 1
state <-  c(B=1278.037,P=64.7532,A=15.04461,H=1713.089,D=100,M=12.75542)
times <- seq(0, 400, by = 1)
out1 <- ode(y= state,times = times, func = CF_meso,parms= parameters, events = list( data =reset_Bac),method='impAdams')

#plot(out1)

Hi <- 1713.089
Di <- 100
LF <- 100*(out1[,'H']/out1[,'D'])/(Hi/Di)
PH <- out1[,'H']*100/Hi

# recriação da figura 5.
par(mfrow=c(2,2))

plot(x=out1[,'time'],y=out0[,'D'],col='black',xlab='time',ylab="",type='l',ylim = c(75,120) )
lines(LF, col='red',lwd=1)
lines(PH, col='blue',lwd=1)
legend('bottomleft', legend=c('D','LF','PH'),
       col=c('black','red','blue'),horiz= TRUE,lty=1,cex = 0.74)

plot(x=out1[,'time'],y=out1[,'B'],col='black',xlab='time',ylab = 'B',type='l')

plot(x=out1[,'time'],y=out1[,'A'],col='black',xlab='time',ylab="",type='l',ylim=c(10,32))
lines(out1[,'M'], col='red',lwd=1)
legend('topleft', legend=c('A','M'),
       col=c('black','red'),lty=1,horiz=TRUE,cex=0.8)

plot(x=out1[,'time'],y=out1[,'P'],col='black',xlab='time',ylab = 'P',type='l')


## Tratamento antibacteriano

# introdução do parametro AB para modelar o tratamento
# AB=5 quando o tratamento está a ser administrado, acontece sempre 2 dias pós infeção, 
# dura 10 dias.

times <- seq(0, 400, by = 1)
signal <- data.frame(times = times, import = rep(0, length(times)))
signal[23:33,2] <- 5
signal[163:173,2] <- 5
signal[303:313,2] <- 5
input <- approxfun(signal, rule = 2)

# modelo CF + tratamento antibacteriano
CF_meso_treat_anti<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    import_AB <- input(t)
    dB = ((B^0.8) * (M^1.2) - (0.8*(B^1.2) * (P^0.1) * (import_AB+1)))*M^- 2.4
    dP <- 20*(B^0.2) * (D^0.4) * (A^(-0.1))-(50*P^0.5)
    dA <- 15*(B^0.1) * (P^0.1)-(12*A^0.5)
    dH <- 500-40*(H^0.3) *  (P^0.2) * (A^(- 0.2))
    dD <- (40*(H^0.3) * (P^0.2) * (A^( -0.2)))-5*D 
    dM = 0.16*(2^CF) * (B^0.25) - 0.15*M^(2-CF)
    # return the rate of change
    list(c(dB,dP,dA,dH,dD,dM))
  }) # end with(as.list ...
}

parameters <- c(CF=1)
state <-  c(B=1278.037,P=64.7532,A=15.04461,H=1713.089,D=100,M=12.75542)
times <- seq(0, 400, by = 1)
out_anti <- ode(y= state,times = times, func = CF_meso_treat_anti ,parms= parameters, events = list(data =reset_Bac),method='impAdams')
#plot(out_anti)

Hi <- 1713.089
Di <- 100
LF <- 100*(out_anti[,'H']/out_anti[,'D'])/(Hi/Di)
PH <- out_anti[,'H']*100/Hi

# recriação da figura 6.
par(mfrow=c(2,2))

plot(x=out_anti[,'time'],y=out_anti[,'D'],col='black',xlab='time',ylab="",type='l',ylim = c(80,120) )
lines(LF, col='red',lwd=1)
lines(PH, col='blue',lwd=1)
legend('bottomleft', legend=c('D','LF','PH'),
       col=c('black','red','blue'),horiz= TRUE,lty=1,cex = 0.74)

plot(x=out_anti[,'time'],y=out_anti[,'B'],col='black',xlab='time',ylab = 'B',type='l')

plot(x=out_anti[,'time'],y=out_anti[,'A'],col='black',xlab='time',ylab="",type='l',ylim=c(10,32))
lines(out_anti[,'M'], col='red',lwd=1)
legend('topleft', legend=c('A','M'),
       col=c('black','red'),lty=1,horiz=TRUE,cex=0.8)


plot(x=out_anti[,'time'],y=out_anti[,'P'],col='black',xlab='time',ylab = 'P',type='l')

# Tratamento antibacteriano + tratamento chest vest 
# introdução do parametro AB para modelar o tratamento chest vest
# Aqui o tratamento CV é considerado como contínuo, começando a t=100
# CV=1 quando o tratamento está 'ativo'

times <- seq(0, 400, by = 1)
signal <- data.frame(times = times, import = rep(0, length(times)))
signal$import <- ifelse(signal$times >100, 1, 0)
input2 <- approxfun(signal, rule = 2)

# modelo CF + tratamento antibacteriano + tratamento CVs
CF_meso_treat_anti_CV<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    import_AB <- input(t)
    import_CV <- input2(t)
    dB = ((B^0.8) * (M^1.2) - (0.8*(B^1.2) * (P^0.1) * (import_AB+1))-0.1*B*import_CV)*M^- 2.4
    dP <- 20*(B^0.2) * (D^0.4) * (A^(-0.1))-(50*P^0.5)
    dA <- 15*(B^0.1) * (P^0.1)-(12*A^0.5)
    dH <- 500-40*(H^0.3) *  (P^0.2) * (A^(- 0.2))
    dD <- (40*(H^0.3) * (P^0.2) * (A^( -0.2)))-5*D 
    dM = 0.16*(2^CF) * (B^0.25) - 0.15*M^(2-CF) - 0.1*M*import_CV
    # return the rate of change
    list(c(dB,dP,dA,dH,dD,dM))
  }) # end with(as.list ...
}

parameters <- c(CF=1)
state <-  c(B=1278.037,P=64.7532,A=15.04461,H=1713.089,D=100,M=12.75542)
times <- seq(0, 400, by = 1)
out_anti_CV <- ode(y= state,times = times, func = CF_meso_treat_anti_CV ,parms= parameters, events = list(data =reset_Bac),method='impAdams')
#plot(out_anti_CV)

Hi <- 1713.089
Di <- 100
LF <- 100*(out_anti_CV[,'H']/out_anti_CV[,'D'])/(Hi/Di)
PH <- out_anti_CV[,'H']*100/Hi

# recriação da fígura 7
par(mfrow=c(2,2))

plot(x=out_anti_CV[,'time'],y=out_anti_CV[,'D'],col='black',xlab='time',ylab="",type='l',ylim = c(80,120) )
lines(LF, col='red',lwd=1)
lines(PH, col='blue',lwd=1)
legend('bottomleft', legend=c('D','LF','PH'),
       col=c('black','red','blue'),horiz= TRUE,lty=1,cex = 0.74)

plot(x=out_anti_CV[,'time'],y=out_anti_CV[,'B'],col='black',xlab='time',ylab = 'B',type='l')

  plot(x=out_anti_CV[,'time'],y=out_anti_CV[,'A'],col='black',xlab='time',ylab="",type='l',ylim=c(0,30))
lines(out_anti_CV[,'M'], col='red',lwd=1)
legend('topleft', legend=c('A','M'),
       col=c('black','red'),lty=1,horiz=TRUE,cex=0.8)


plot(x=out_anti_CV[,'time'],y=out_anti_CV[,'P'],col='black',xlab='time',ylab = 'P',type='l')


# Expansão da variável B para uma metapopulação de 3 espécies de bacterias
### População de 3 espécies de bacterias 

# modelo de interação 3 especies bacterianas 
# B aqui é uma aproximado (simplificação) da variável B do modelo da doença CF
CF_bac <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dB1 <- 0.01*B1 - 1.2*10^-5*B1^2 - 6*10^-6 * B1 * B3
    dB2 <- 0.03 * B2^0.8 + 1.75 * 10^-5 * B1 * B2 - 0.002*B2^1.2 - 2.0 * 10^-5 * B2 * B3
    dB3 <- 0.06*B3^0.5 - 2.0 * 10^-9 * B3^3
    dB <- 0.047*B^0.8-0.0027*B^1.2
    # return the rate of change
    list(c(dB,dB1,dB2,dB3))
  }) # end with(as.list ...
}


# valores iniciais de B,B1,B2 e B3
state <-  c(B=1,B1=1,B2=1.2,B3=1)
times <- seq(0, 2500, by = 1)
out_bac <- ode(y= state,times = times, func = CF_bac ,parms= NULL,method='impAdams')
#plot(out_bac)

# soma de B1+B2+B2
library(dplyr)
out_bac_df <- as.data.frame(out_bac)
out_bac_df <- mutate(out_bac_df,sum=B1+B2+B3)

# recriação da fígura 9
# legenda : azul=S, preto=B, roxo=B3, verde=B2, vermelho=B1.
# (a legenda ficava ilegível)

par(mfrow=c(2,2))

# plot 9a: crescimento das espécies até ao steady state + soma das três 
# e comparação com o valor de B
plot(x=out_bac_df$time,y=out_bac_df$sum,col='blue',type='l',xlab='time',ylab = 'size',main = 'bac pop',lwd=1)
lines(out_bac_df$B3, col='purple',lwd=1)
lines(out_bac_df$B2, col='red',lwd=1)
lines(out_bac_df$B1, col='green',lwd=1)
lines(out_bac_df$B, col='black',lwd=1)
#legend('topleft', legend=c('sum','B3','B2','B1','B'),
 #      col=c('green','yellow','blue','red','black'),lty=1)


# df de eventos: matamos 99% da população de cada espécie
# quando t=3000
kill99_Bac <- data.frame( 
  var = c('B','B1','B2','B3'), 
  time = 3000 , 
  value = 0.01,
  method = 'multiply')

state <-  c(B=1,B1=1,B2=1.2,B3=1)
times <- seq(0, 5000, by = 1)
out_kill99_bac <- ode(y= state,times = times, func = CF_bac,parms= NULL, events = list(data=kill99_Bac),method='impAdams')
#plot(out_kill99_bac)

out_kill99_bac_df <- as.data.frame(out_kill99_bac)
out_kill99_bac_df <- mutate(out_kill99_bac_df,sum=B1+B2+B3)

# plot 9b : matamos 99% da população de cada espécie
plot(x=out_kill99_bac_df$time,y=out_kill99_bac_df$sum,type='l',col='blue',xlab='time',ylab = 'size',main = 'bac pop',lwd=1,ylim=c(0,2000))
lines(out_kill99_bac_df$B3, col='purple',lwd=1)
lines(out_kill99_bac_df$B2, col='red',lwd=1)
lines(out_kill99_bac_df$B1, col='green',lwd=1)
lines(out_kill99_bac_df$B, col='black',lwd=1)
#legend('topleft', legend=c('sum','B3','B2','B1','B'),
 #      col=c('green','purple','blue','red','black',horiz=TRUE),lty=1)

## df com eventos. Matamos 99% da população de B2 e B3 e 100% de B1.
killB1_Bac <- data.frame( 
  var = c('B1','B2','B3'), 
  time = 5500 , 
  value = c(0,0.01,0.01),
  method = 'multiply')

state <-  c(B=1,B1=1,B2=1.2,B3=1)
times <- seq(0, 10000, by = 1)
out_killB1_bac <- ode(y= state,times = times, func = CF_bac,parms= NULL, events = list(data=killB1_Bac),method='impAdams')
#plot(out_killB1_bac)

# plot 9c : Matamos 99% da população de B2 e B3 e 100% de B1.
plot(x=out_killB1_bac[,'time'],y=out_killB1_bac[,'B1'],type='l',col='green',xlab='time',ylab = 'size',main = 'bac pop',lwd=1,ylim=c(0,1000))
lines(out_killB1_bac[,'B2'], col='red',lwd=1)
lines(out_killB1_bac[,'B3'], col='purple',lwd=1)
#legend('topleft', legend=c('B1','B2','B3'),
 #      col=c('green','red','blue'),lty=1)

## df com eventos. Matamos 99% da população de B1 e B2 e 100% de B3.
killB3_Bac <- data.frame( 
  var = c('B1','B2','B3'), 
  time = 5500 , 
  value = c(0.01,0.01,0),
  method = 'multiply')

state <-  c(B=1,B1=1,B2=1.2,B3=1)
times <- seq(0, 10000, by = 1)
out_killB3_bac <- ode(y= state,times = times, func = CF_bac,parms= NULL, events = list(data=killB3_Bac),method='impAdams')
#plot(out_killB3_bac)

# plot 9d : Matamos 99% da população de B1 e B2 e 100% de B3.
plot(x=out_killB3_bac[,'time'],y=out_killB3_bac[,'B1'],type='l',col='green',xlab='time',ylab = 'size',main = 'bac pop',lwd=1,ylim=c(0,60000))
lines(out_killB3_bac[,'B2'], col='red',lwd=1)
lines(out_killB3_bac[,'B3'], col='purple',lwd=1)
#legend('topleft', legend=c('B1','B2','B3'),
 #      col=c('green','red','blue'),lty=1)


##
### Não mostrado mas mencionado :

# matar 100% de B2 e 99% de B1 e B3
killB2_Bac <- data.frame( 
  var = c('B1','B2','B3'), 
  time = 5500 , 
  value = c(0.01,0,0.01),
  method = 'multiply')

plot.new()
state <-  c(B=1,B1=1,B2=1.2,B3=1)
times <- seq(0, 10000, by = 1)
out_killB2_bac <- ode(y= state,times = times, func = CF_bac,parms= NULL, events = list(data=killB2_Bac),method='impAdams')
#Matamos 99% da população de B1 e B3 e 100% de B2.
plot(x=out_killB2_bac[,'time'],y=out_killB2_bac[,'B1'],type='l',col='green',xlab='time',ylab = 'size',main = 'bac pop',lwd=1,ylim=c(0,1000))
lines(out_killB2_bac[,'B2'], col='red',lwd=1)
lines(out_killB2_bac[,'B3'], col='purple',lwd=1)

# 

# Diferentes interações de bacterias levam a diferentes
# dinamicas de sistemas, como por exemplo osciladores:

CF_bac_osc <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dB1 <- 0.035*(B2-B1^0.1)
    dB2 <- 0.035*((B1^-0.5)*B2^0.2-B2^0.1)
    dB3 <- 0.8*((B2^-0.2)*B3 - B1*B3^2)
    list(c(dB1,dB2,dB3))
  }) # end with(as.list ...
}

state <-  c(B1=1,B2=1.2,B3=1)
times <- seq(0, 1000, by = 1)
out_bac_osc <- ode(y= state,times = times, func = CF_bac_osc,parms= NULL,method='impAdams')
#plot(out_bac_osc)

# recriação da fígura 10
plot.new()
par(mfrow=c(1,1))

plot(x=out_bac_osc[,'time'],y=out_bac_osc[,'B1'],type='l',col='green',xlab='time',ylab = 'size',main = 'bac pop',lwd=1,ylim=c(0,2))
lines(out_bac_osc[,'B2'], col='red',lwd=1)
lines(out_bac_osc[,'B3'], col='purple',lwd=1)
legend('topleft', legend=c('B1','B2','B3'),
       col=c('green','red','purple'),lty=1)

# Expansão das variáveis P e A em conjunto. 
### Modelo interação P e A

CF_PA <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dP <- 0.1*(1+miu)*((A^-0.15)*P^0.4-P^0.2)
    dA <- 0.1*(P-A^0.2)
    list(c(dP,dA))
  }) # end with(as.list ...
}

# oscilações 'damped' miu = -0.2
parameters <- c(miu=-0.2)
# valores iniciais de P e A
state <-  c(P=1,A=1.8)
times <- seq(0, 3000, by = 1)
out_PA_damped <- ode(y= state,times = times, func = CF_PA,parms= parameters,method='impAdams')

#plot(out_PA_damped)

# recriação da fígura 12
par(mfrow=c(1,2))

# oscilações 'damped'
plot(x=out_PA_damped[,'time'],y=out_PA_damped[,'P'],type='l',col='blue',xlab='time',ylab = '',main = 'PA',lwd=1,ylim=c(0,2))
lines(out_PA_damped[,'A'], col='red',lwd=1)
legend('topleft', legend=c('P','A'),
       col=c('blue','red'),lty=1)



# oscilações 'stable limit cyle' miu = 0.04
# resultado semelhante se miu = 0.001
parameters <- c(miu=0.04)
times <- seq(0, 1000, by = 1)
out_PA_osc <- ode(y= state,times = times, func = CF_PA,parms= parameters,method='impAdams')
#plot(out_PA_osc)

plot(x=out_PA_osc[,'time'],y=out_PA_osc[,'P'],type='l',col='blue',xlab='time',ylab = '',main = 'PA',lwd=1,ylim=c(0,2))
lines(out_PA_osc[,'A'], col='red',lwd=1)
#legend('topleft', legend=c('P','A'),
 #      col=c('blue','red'),lty=1)


### modelo metapoplação bacteriano + P e A

CF_PA_bac_osc <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dB1 <- 0.035*(B2-B1^0.1)
    dB2 <- 0.035*((B1^-0.5)*B2^0.2-B2^0.1)
    dB3 <- 0.8*((B2^-0.2)*B3 - B1*B3^2)
    dP <- 0.1*(1+0.04)*((A^-0.15)*P^0.4+q*B3-P^0.2-q)
    dA <- 0.1*(P-A^0.2)
    list(c(dB1,dB2,dB3,dP,dA))
  }) # end with(as.list ...
}

# oscilação irregular 
parameters <- c(q=0.01)
state <-  c(B1=1,B2=1.2,B3=1,P=1,A=1.8)
times <- seq(0, 1000, by = 1)
out__bac_PA_001 <- ode(y= state,times = times, func = CF_PA_bac_osc,parms= parameters,method='impAdams')
#plot(out__bac_PA_001)

# recriação da fígura 13
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

plot(x=out__bac_PA_001[,'time'],y=out__bac_PA_001[,'P'],type='l',col='blue',xlab='time',ylab = '',main = 'PA',lwd=1,ylim=c(0,2))
lines(out__bac_PA_001[,'A'], col='red',lwd=1)
#legend('topright', legend=c('P','A'),
 #      col=c('blue','red'),lty=1,horiz=TRUE)

# oscilação dupla irregular
parameters <- c(q=0.2)
out__bac_PA_02 <- ode(y= state,times = times, func = CF_PA_bac_osc,parms= parameters,method='impAdams')
#plot(out__bac_PA_02)

plot(x=out__bac_PA_02[,'time'],y=out__bac_PA_02[,'P'],type='l',col='blue',xlab='time',ylab = '',main = 'PA',lwd=1,ylim=c(0,2.2))
lines(out__bac_PA_02[,'A'], col='red',lwd=1)
#legend('topleft', legend=c('P','A'),
 #      col=c('blue','red'),lty=1)

# dinâmica errática, possivelmente caótica.
parameters <- c(q=0.08)
times <- seq(0, 20000, by = 1)
out__bac_PA_008 <- ode(y= state,times = times, func = CF_PA_bac_osc,parms= parameters,method='impAdams')

plot(x=out__bac_PA_008[,'time'],y=out__bac_PA_008[,'A'],type='l',col='blue',xlab='time',ylab = '',main = 'PA',lwd=1,xlim=c(15000,20000))
lines(out__bac_PA_008[,'P'], col='red',lwd=1)
#legend('topleft', legend=c('A','P'),
 #      col=c('blue','red'),lty=1)

# recriação da fígura 14. P x A
par(mfrow=c(1,1))
plot(x=out__bac_PA_008[,'A'],y=out__bac_PA_008[,'P'],type='l',col='blue',xlab='A',ylab = 'P',main = 'PA',lwd=1)

