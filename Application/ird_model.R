library(magidde)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1)
### use 30-day data in parameter/trajectory inference ####
sirmodel = list(name = "IRD",
                 thetaLowerBound = c(0, 0, 0, 0),
                 thetaUpperBound = c(Inf, Inf, Inf, Inf))
N = 14999441
dat = read.csv("sir_data.csv", header = TRUE)[1:30,2:5]
y = setDiscretization(dat, level  = 1)
res1_sirtwo = MagiSolver(y, sirmodel, control = list(nstepsHmc = 25, 
                                                     niterHmc = 80000))
knitr::kable(summary(res1_sirtwo, sigma = TRUE))

### Generate Figure ####
tvec = res1_sirtwo$tvec
t_obs = seq(0, 29, 1)
tlen = length(res1_sirtwo$tvec)


df_I = data.frame(time = tvec,
                 component = rep("I", tlen),
                 est = colMeans(res1_sirtwo$xsampled[,,1]),
                 q.025 =  apply(res1_sirtwo$xsampled[,,1], MARGIN = 2, function (x) quantile(x, 0.025)),
                 q.975 =  apply(res1_sirtwo$xsampled[,,1], MARGIN = 2, function (x) quantile(x, 0.975)))

df_R = data.frame(time = tvec,
                 component = rep("R", tlen),
                 est = colMeans(res1_sirtwo$xsampled[,,2]),
                 q.025 =  apply(res1_sirtwo$xsampled[,,2], MARGIN = 2, function (x) quantile(x, 0.025)),
                 q.975 =  apply(res1_sirtwo$xsampled[,,2], MARGIN = 2, function (x) quantile(x, 0.975)))

df_D = data.frame(time = tvec,
                 component = rep("D", tlen),
                 est = colMeans(res1_sirtwo$xsampled[,,3]),
                 q.025 =  apply(res1_sirtwo$xsampled[,,3], MARGIN = 2, function (x) quantile(x, 0.025)),
                 q.975 =  apply(res1_sirtwo$xsampled[,,3], MARGIN = 2, function (x) quantile(x, 0.975)))


df_S = data.frame(time = tvec,
                 component = rep("S", tlen),
                 est = 1-df_I$est -df_R$est-df_D$est,
                 q.025 = 1-df_I$q.025 -df_R$q.025-df_D$q.025,
                 q.975 = 1-df_I$q.975 -df_R$q.975-df_D$q.975)



df_I$est = N * df_I$est
df_I$q.025 = N * df_I$q.025
df_I$q.975 = N * df_I$q.975


df_R$est = N * df_R$est
df_R$q.025 = N * df_R$q.025
df_R$q.975 = N * df_R$q.975


df_D$est = N * df_D$est
df_D$q.025 = N * df_D$q.025
df_D$q.975 = N * df_D$q.975

df_S$est = N * df_S$est
df_S$q.025 = N * df_S$q.025
df_S$q.975 = N * df_S$q.975


data_I = c(na.omit(res1_sirtwo$y[,1])) * N
data_R = c(na.omit(res1_sirtwo$y[,2])) * N
data_D = c(na.omit(res1_sirtwo$y[,3])) * N
data_S = c(1-(na.omit(res1_sirtwo$y[,3])+na.omit(res1_sirtwo$y[,2])+na.omit(res1_sirtwo$y[,1]))) * N


obs_df = data.frame(time=t_obs,
                    S = data_S,
                    I = data_I,
                    R  = data_R,
                    D = data_D)


p_1 = ggplot()+
  geom_ribbon(data = df_S, aes(ymin = q.025, ymax = q.975, x = time,fill = "95% credible interval"), alpha = 0.4)+
  geom_line(data = df_S,aes(x = time, y = est, colour = "Inferred trajectories"))+
  theme_bw()+ 
  geom_point(data = obs_df, aes(x = t_obs, y = S, shape = "Observations"), colour = "black",size =0.5)+
  scale_fill_manual(values=c("95% credible interval" = "#66c2a4"),labels=c("95% credible interval"),
                    name = "")+
  scale_colour_manual(values = c("Inferred trajectories" = "#006d2c"),name = "")+
  scale_shape_manual(values = c("Observations" = 1),name = "")+ 
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text =  element_text(size = 12), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+xlab("Time")+ylab("S") + 
  guides(fill = guide_legend(order = 2), linetype = guide_legend(order = 1), colour = guide_legend(order = 1))


p_2 = ggplot()+
  geom_ribbon(data = df_I, aes(ymin = q.025, ymax = q.975, x = time,fill = "95% credible interval"), alpha = 0.4)+
  geom_line(data = df_I,aes(x = time, y = est, colour = "Inferred trajectories"))+
  theme_bw()+ 
  geom_point(data = obs_df, aes(x = t_obs, y = I, shape = "Observations"), colour = "black",size =0.5)+
  scale_fill_manual(values=c("95% credible interval" = "#66c2a4"),labels=c("95% credible interval"),
                    name = "")+
  scale_colour_manual(values = c("Inferred trajectories" = "#006d2c"),name = "")+
  scale_shape_manual(values = c("Observations" = 1),name = "")+ 
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text =  element_text(size = 12), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+xlab("Time")+ylab("I") + 
  guides(fill = guide_legend(order = 2), linetype = guide_legend(order = 1), colour = guide_legend(order = 1))


p_3 = ggplot()+
  geom_ribbon(data = df_R, aes(ymin = q.025, ymax = q.975, x = time,fill = "95% credible interval"), alpha = 0.4)+
  geom_line(data = df_R, aes(x = time, y = est, colour = "Inferred trajectories"))+
  theme_bw()+ 
  geom_point(data = obs_df,aes(x = time, y = R, shape = "Observations"), colour = "black",size =0.5  )+
  scale_fill_manual(values=c("95% credible interval" = "#66c2a4"),labels=c("95% credible interval"),
                    name = "")+
  scale_colour_manual(values = c("Inferred trajectories" = "#006d2c"),name = "")+
  scale_shape_manual(values = c("Observations" = 1),name = "")+ 
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text =  element_text(size = 12), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+xlab("Time")+ylab("R") + 
  guides(fill = guide_legend(order = 2), linetype = guide_legend(order = 1), colour = guide_legend(order = 1))

p_4 = ggplot()+
  geom_ribbon(data = df_D,aes(ymin = q.025, ymax = q.975, x = time,fill = "95% credible interval"), alpha = 0.4)+
  geom_line(data = df_D,aes(x = time, y = est, colour = "Inferred trajectories"))+
  theme_bw()+ 
  geom_point(data = obs_df, aes(x = time, y = data_D, shape = "Observations"), colour = "black",size =0.5  )+
  scale_fill_manual(values=c("95% credible interval" = "#66c2a4"),labels=c("95% credible interval"),
                    name = "")+
  scale_colour_manual(values = c("Inferred trajectories" = "#006d2c"),name = "")+
  scale_shape_manual(values = c("Observations" = 1),name = "")+ 
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text =  element_text(size = 12), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+xlab("Time")+ylab("D") + 
  guides(fill = guide_legend(order = 2), linetype = guide_legend(order = 1), colour = guide_legend(order = 1))


p = p_1+p_2+p_3+p_4 + plot_layout(guides = "collect",ncol = 4) & theme(legend.position = 'bottom')
p

### use 16-day data in parameter/trajectory inference ####

dat_train =  read.csv("sir_data.csv", header = TRUE)[,2:5]
dat_train[c(17:31),c(2, 3, 4)] = NA
y_new = setDiscretization(dat_train, level  = 1)

res1_sir_train = MagiSolver(y_new, sirmodel,  
                                     control = list(nstepsHmc = 25, niterHmc = 80000))

### Generate Figure ####

tvec = res1_sir_train$tvec

tlen = length(res1_sir_train$tvec)
train_day = 15



df_I = data.frame(time = tvec,
                  component = rep("I", tlen),
                  est = colMeans(res1_sir_train$xsampled[,,1]),
                  q.025 =  apply(res1_sir_train$xsampled[,,1], MARGIN = 2, function (x) quantile(x, 0.025)),
                  q.975 =  apply(res1_sir_train$xsampled[,,1], MARGIN = 2, function (x) quantile(x, 0.975)),
                  is_train = c(rep("train", train_day * 2 + 1), rep("test", 60 - train_day * 2)))


df_R =data.frame(time = tvec,
                 component = rep("R", tlen),
                 est = colMeans(res1_sir_train$xsampled[,,2]),
                 q.025 =  apply(res1_sir_train$xsampled[,,2], MARGIN = 2, function (x) quantile(x, 0.025)),
                 q.975 =  apply(res1_sir_train$xsampled[,,2], MARGIN = 2, function (x) quantile(x, 0.975)),
                 is_train = c(rep("train", train_day * 2 + 1), rep("test", 60 - train_day * 2)))

df_D =data.frame(time = tvec,
                 component = rep("D", tlen),
                 est = colMeans(res1_sir_train$xsampled[,,3]),
                 q.025 =  apply(res1_sir_train$xsampled[,,3], MARGIN = 2, function (x) quantile(x, 0.025)),
                 q.975 =  apply(res1_sir_train$xsampled[,,3], MARGIN = 2, function (x) quantile(x, 0.975)),
                 is_train = c(rep("train", train_day * 2 + 1), rep("test", 60 - train_day * 2)))




df_S =data.frame(time = tvec,
                 component = rep("S", tlen),
                 est = 1-df_I$est -df_R$est-df_D$est,
                 q.025 = 1-df_I$q.025 -df_R$q.025-df_D$q.025,
                 q.975 = 1-df_I$q.975 -df_R$q.975-df_D$q.975,
                 is_train = c(rep("train", train_day * 2 + 1), rep("test", 60 - train_day * 2)))


df_I <- add_row(df_I, df_I[train_day * 2 + 1,], .after = train_day * 2 + 1)
df_I[train_day * 2 + 2, "is_train"] = "test"
df_R <- add_row(df_R, df_R[train_day * 2 + 1,], .after = train_day * 2 + 1)
df_R[train_day * 2 + 2, "is_train"] = "test"
df_D <- add_row(df_D, df_D[train_day * 2 + 1,], .after = train_day * 2 + 1)
df_D[train_day * 2 + 2, "is_train"] = "test"
df_S <- add_row(df_S, df_S[train_day * 2 + 1,], .after = train_day * 2 + 1)
df_S[train_day * 2 + 2, "is_train"] = "test"


df_I$est = N * df_I$est
df_I$q.025 = N * df_I$q.025
df_I$q.975 = N * df_I$q.975


df_R$est = N * df_R$est
df_R$q.025 = N * df_R$q.025
df_R$q.975 = N * df_R$q.975


df_D$est = N * df_D$est
df_D$q.025 = N * df_D$q.025
df_D$q.975 = N * df_D$q.975

df_S$est = N * df_S$est
df_S$q.025 = N * df_S$q.025
df_S$q.975 = N * df_S$q.975


sir_data = read.csv("sir_data.csv")
data_df = data.frame(time=seq(0, 30, 1),
                     S =(1 - (sir_data$I[1:31] + sir_data$R[1:31] + sir_data$D[1:31]))* N,
                     I = sir_data$I[1:31]*N,
                     R  = sir_data$R[1:31]*N,
                     D = sir_data$D[1:31]*N,
                     is_train = c(rep("train", train_day+1), rep("test", 30-train_day)))


train_fill <- "#66c2a4" 
predict_fill <- "#fdae6b" 
line_color <- "#006d2c" 


df_I$is_train <- factor(df_I$is_train, levels = c("train", "test"))
data_df$is_train <- factor(data_df$is_train, levels = c("train", "test"))
df_S$is_train <- factor(df_S$is_train, levels = c("train", "test"))
df_R$is_train <- factor(df_R$is_train, levels = c("train", "test"))
df_D$is_train <- factor(df_D$is_train, levels = c("train", "test"))

df_S

p_1 <- ggplot() +
  geom_ribbon(data = df_S, aes(ymin = q.025, ymax = q.975, x = time, group = is_train, fill = is_train), alpha = 0.4) +
  geom_line(data = df_S, aes(x = time, y = est, color = "Inferred trajectories")) +
  geom_point(data = data_df, aes(x = time, y = S, shape = is_train), color = "black", size = 3) +
  scale_color_manual(values = c("Inferred trajectories" = line_color),
                     labels = c("Inferred trajectories" = "Inferred trajectories"))+
  scale_fill_manual(values = c("train" = train_fill, "test" = predict_fill), 
                    labels = c("train" = "95% credible interval in fitting", "test" = "95% credible interval in prediction"), name = "95% Credible Interval") +
  scale_shape_manual(values = c("train" = 16, "test" = 1), 
                     labels = c("train" = "Observations for fitting", "test" = "Observations for prediction"), name = "Observations") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab(" ") + ylab("S") + geom_vline(xintercept = train_day, linetype = "dashed")



p_2 <- ggplot() +
  geom_ribbon(data = df_I, aes(ymin = q.025, ymax = q.975, x = time, group = is_train, fill = is_train), alpha = 0.4) +
  geom_line(data = df_I, aes(x = time, y = est, color = "Inferred trajectories")) +
  geom_point(data = data_df, aes(x = time, y = I, shape = is_train), color = "black", size = 3) +
  scale_color_manual(values = c("Inferred trajectories" = line_color),
                     labels = c("Inferred trajectories" = "Inferred trajectories"))+
  scale_fill_manual(values = c("train" = train_fill, "test" = predict_fill), 
                    labels = c("train" = "95% credible interval in fitting", "test" = "95% credible interval in prediction"), name = "95% Credible Interval") +
  scale_shape_manual(values = c("train" = 16, "test" = 1), 
                     labels = c("train" = "Observations for fitting", "test" = "Observations for prediction"), name = "Observations") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab(" ") + ylab("I") + geom_vline(xintercept = train_day, linetype = "dashed")


p_3 <- ggplot() +
  geom_ribbon(data = df_R, aes(ymin = q.025, ymax = q.975, x = time, group = is_train, fill = is_train), alpha = 0.4) +
  geom_line(data = df_R, aes(x = time, y = est, color = "Inferred trajectories")) +
  geom_point(data = data_df, aes(x = time, y = R, shape = is_train), color = "black", size = 3) +
  scale_color_manual(values = c("Inferred trajectories" = line_color),
                     labels = c("Inferred trajectories" = "Inferred trajectories"))+
  scale_fill_manual(values = c("train" = train_fill, "test" = predict_fill), 
                    labels = c("train" = "95% credible interval in fitting", "test" = "95% credible interval in prediction"), name = "95% Credible Interval") +
  scale_shape_manual(values = c("train" = 16, "test" = 1), 
                     labels = c("train" = "Observations for fitting", "test" = "Observations for prediction"), name = "Observations") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab("Time") + ylab("R") + geom_vline(xintercept = train_day, linetype = "dashed")

p_4 <- ggplot() +
  geom_ribbon(data = df_D, aes(ymin = q.025, ymax = q.975, x = time, group = is_train, fill = is_train), alpha = 0.4) +
  geom_line(data = df_D, aes(x = time, y = est, color = "Inferred trajectories")) +
  geom_point(data = data_df, aes(x = time, y = D, shape = is_train), color = "black", size = 3) +
  scale_color_manual(values = c("Inferred trajectories" = line_color),
                     labels = c("Inferred trajectories" = "Inferred trajectories"))+
  scale_fill_manual(values = c("train" = train_fill, "test" = predict_fill), 
                    labels = c("train" = "95% credible interval in fitting", "test" = "95% credible interval in prediction"), name = "95% Credible Interval") +
  scale_shape_manual(values = c("train" = 16, "test" = 1), 
                     labels = c("train" = "Observations for fitting", "test" = "Observations for prediction"), name = "Observations") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab("Time") + ylab("D") + geom_vline(xintercept = train_day, linetype = "dashed")

# Combine plots
p <- p_1 + p_2 + p_3 + p_4 + plot_layout(guides = "collect", ncol = 2) & 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.box = "horizontal") 

p

