library(deSolve)
library(magidde)
library(patchwork)
library(dplyr)
library(ggplot2)
tvec = seq(0, 30, 1)
N = 14999441
tlen = length(tvec)
set.seed(1)
sir_dde = function(t, y, parms){
  beta = parms[1]
  h = parms[2]
  mu_d = parms[3]
  lambda = parms[4]
  I_init = parms[5]
  
  S = y[1]
  I = y[2]
  R = y[3]
  D = y[4]
  
  if (t < h) lag_I = I_init 
  else lag_I = lagvalue(t-h)[2]  
  dS = -beta * S * lag_I
  dI = beta * S * lag_I - mu_d * I - lambda * I
  dR = lambda * I
  dD = mu_d * I 
  
  return(list(c(dS, dI, dR, dD)))
}

y_full = dede(func = sir_dde, times = tvec, parms = c(0.02542701, 3.035970, 0.0003326903, 
                                                      0.07508062, 0.014510457), 
              y = c(1-(0.014510457 + 0.005306228 + 7.484958e-05),0.014510457, 0.005306228,  7.484958e-05))
y_full_IRD = data.frame(y_full[,c(1, 3, 4, 5)])
sigma_vec = c(3.473044e-04, 3.475568e-04,  2.155736e-07)
y_obs = data.frame(cbind(y_full_IRD$time, c(y_full_IRD$X2 + rnorm(tlen, 0, sigma_vec[1])),
                         c(y_full_IRD$X3 + rnorm(tlen, 0, sigma_vec[2])),
                         c(y_full_IRD$X4 + rnorm(tlen, 0, sigma_vec[3]))))
colnames(y_obs)  = c("time", "I", "R", "D")
y_obs$S = 1 - (y_obs$I + y_obs$R + y_obs$D)
y_obs_train = y_obs[,c(1:4)]
y_obs_train[17:31, 2:4] = NA 
y_In = setDiscretization(y_obs_train, level = 1)
sirmodel = list(name = "IRD",
                thetaLowerBound = c(0, 0, 0, 0),
                thetaUpperBound = c(Inf, Inf, Inf, Inf))
res1_sir = MagiSolver(y_In, sirmodel, control = list(nstepsHmc = 25,
                                                     niterHmc = 80000))

summary(res1_sir, sigma = TRUE)

### Generate figure ####

train_day = 15
tvec_discretization =  res1_sir$tvec
tlen_discretization = length(tvec_discretization)
df_I = data.frame(time = tvec_discretization,
                  component = rep("I", tlen_discretization),
                  est = colMeans(res1_sir$xsampled[,,1]),
                  q.025 =  apply(res1_sir$xsampled[,,1], MARGIN = 2, function (x) quantile(x, 0.025)),
                  q.975 =  apply(res1_sir$xsampled[,,1], MARGIN = 2, function (x) quantile(x, 0.975)),
                  is_train = c(rep("train", train_day * 2 + 1), rep("test", 60 - train_day * 2)))

df_R = data.frame(time = tvec_discretization,
                  component = rep("R", tlen_discretization),
                  est = colMeans(res1_sir$xsampled[,,2]),
                  q.025 =  apply(res1_sir$xsampled[,,2], MARGIN = 2, function (x) quantile(x, 0.025)),
                  q.975 =  apply(res1_sir$xsampled[,,2], MARGIN = 2, function (x) quantile(x, 0.975)),
                  is_train = c(rep("train", train_day * 2 + 1), rep("test", 60 - train_day * 2)))

df_D = data.frame(time = tvec_discretization,
                  component = rep("D", tlen_discretization),
                  est = colMeans(res1_sir$xsampled[,,3]),
                  q.025 =  apply(res1_sir$xsampled[,,3], MARGIN = 2, function (x) quantile(x, 0.025)),
                  q.975 =  apply(res1_sir$xsampled[,,3], MARGIN = 2, function (x) quantile(x, 0.975)),
                  is_train = c(rep("train", train_day * 2 + 1), rep("test", 60 - train_day * 2)))



df_S =data.frame(time = tvec_discretization,
                 component = rep("S", tlen_discretization),
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

data_df_truth= data.frame(time = tvec,
                          S = y_full[,2] * N,
                          I = y_full[,3] * N,
                          R  = y_full[,4]* N,
                          D = y_full[,5] * N,
                          is_train = c(rep("train", train_day+1), rep("test", 30-train_day)))


data_df_obs = data.frame(time=seq(0, 30, 1),
                         S = y_obs[,5]* N,
                         I = y_obs[,2]*N,
                         R  = y_obs[,3]*N,
                         D = y_obs[,4]*N,
                         is_train = c(rep("train", train_day+1), rep("test", 30-train_day)))


train_fill <- "#66c2a4" 
predict_fill <- "#fdae6b" 
line_color <- "#006d2c" 


df_I$is_train <- factor(df_I$is_train, levels = c("train", "test"))
data_df_truth$is_train <- factor(data_df_truth$is_train, levels = c("train", "test"))
data_df_obs$is_train <- factor(data_df_obs$is_train, levels = c("train", "test"))
df_S$is_train <- factor(df_S$is_train, levels = c("train", "test"))
df_R$is_train <- factor(df_R$is_train, levels = c("train", "test"))
df_D$is_train <- factor(df_D$is_train, levels = c("train", "test"))


df_S$is_train_S <- df_S$is_train
data_df_truth$is_train_truth <- data_df_truth$is_train
levels(data_df_truth$is_train_truth) <- c("train_truth", "test_truth")
data_df_obs$is_train_obs <- data_df_obs$is_train
levels(data_df_obs$is_train_obs) <- c("train_obs", "test_obs")


df_I$is_train_I <- df_I$is_train
df_R$is_train_R <- df_R$is_train
df_D$is_train_D <- df_D$is_train

p_1 <- ggplot() +
  geom_ribbon(data = df_S, aes(ymin = q.025, ymax = q.975, x = time, fill = is_train_S), alpha = 0.4) +
  geom_line(data = data_df_truth, aes(x = time, y = S, color = "True trajectories"), size = 3) +
  geom_line(data = df_S, aes(x = time, y = est, color = "Inferred trajectories"), size = 1.75) +
  geom_point(data = data_df_obs, aes(x = time, y = S, shape = is_train_obs, color = is_train_obs), size = 3) +
  scale_fill_manual(values = c("train" = train_fill, "test" = predict_fill), 
                    labels = c("train" = "95% credible interval in fitting", 
                               "test" = "95% credible interval in prediction"), 
                    name = NULL) +
  scale_shape_manual(values = c("Inferred trajectories" = NA, 
                                "train_obs" = 16, "test_obs" = 1), 
                     labels = c("Inferred trajectories" = "Inferred trajectories",
                                "train_obs" = "Observations for fitting", 
                                "test_obs" = "Observations for prediction"),
                     name = NULL) +
  scale_color_manual(values = c("Inferred trajectories" = line_color,
                                "True trajectories" = "red",
                                "train_obs" = "black", "test_obs" = "black"), 
                     labels = c("Inferred trajectories" = "Inferred trajectories",
                                "True trajectories" = "True trajectories",
                                "train_obs" = "Observations for fitting", 
                                "test_obs" = "Observations for prediction"
                     ),
                     name = NULL) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape = FALSE,
         col =guide_legend(override.aes = list(linetype = c(0,0,1,1),
                                               shape = c(16,1,NA,NA),
                                               color = c("black" , "black",
                                                         line_color,
                                                         "red"))))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab(NULL) + ylab("S") + geom_vline(xintercept = train_day, linetype = "dashed")



p_2 <- ggplot() +
  geom_ribbon(data = df_I, aes(ymin = q.025, ymax = q.975, x = time, fill = is_train_I), alpha = 0.4) +
  geom_line(data = data_df_truth, aes(x = time, y = I, color = "True trajectories"), size = 3) +
  geom_line(data = df_I, aes(x = time, y = est, color = "Inferred trajectories"), size = 1.75) +
  geom_point(data = data_df_obs, aes(x = time, y = I, shape = is_train_obs, color = is_train_obs), size = 3) +
  scale_fill_manual(values = c("train" = train_fill, "test" = predict_fill), 
                    labels = c("train" = "95% credible interval in fitting", 
                               "test" = "95% credible interval in prediction"), 
                    name = NULL) +
  scale_shape_manual(values = c("Inferred trajectories" = NA, 
                                "train_obs" = 16, "test_obs" = 1), 
                     labels = c("Inferred trajectories" = "Inferred trajectories",
                                "train_obs" = "Observations for fitting", 
                                "test_obs" = "Observations for prediction"),
                     name = NULL) +
  scale_color_manual(values = c("Inferred trajectories" = line_color,
                                "True trajectories" = "red",
                                "train_obs" = "black", "test_obs" = "black"), 
                     labels = c("Inferred trajectories" = "Inferred trajectories",
                                "True trajectories" = "True trajectories",
                                "train_obs" = "Observations for fitting", 
                                "test_obs" = "Observations for prediction"
                     ),
                     name = NULL) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape = FALSE,
         col =guide_legend(override.aes = list(linetype = c(0,0,1,1),
                                               shape = c(16,1,NA,NA),
                                               color = c("black" , "black",
                                                         line_color,
                                                         "red"))))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab(NULL) + ylab("I") + geom_vline(xintercept = train_day, linetype = "dashed")



p_3 <- ggplot() +
  geom_ribbon(data = df_R, aes(ymin = q.025, ymax = q.975, x = time, fill = is_train_R), alpha = 0.4) +
  geom_line(data = data_df_truth, aes(x = time, y = R, color = "True trajectories"), size = 3) +
  geom_line(data = df_R, aes(x = time, y = est, color = "Inferred trajectories"), size = 1.75) +
  geom_point(data = data_df_obs, aes(x = time, y = R, shape = is_train_obs, color = is_train_obs), size = 3) +
  scale_fill_manual(values = c("train" = train_fill, "test" = predict_fill), 
                    labels = c("train" = "95% credible interval in fitting", 
                               "test" = "95% credible interval in prediction"), 
                    name = NULL) +
  scale_shape_manual(values = c("Inferred trajectories" = NA, 
                                "train_obs" = 16, "test_obs" = 1), 
                     labels = c("Inferred trajectories" = "Inferred trajectories",
                                "train_obs" = "Observations for fitting", 
                                "test_obs" = "Observations for prediction"),
                     name = NULL) +
  scale_color_manual(values = c("Inferred trajectories" = line_color,
                                "True trajectories" = "red",
                                "train_obs" = "black", "test_obs" = "black"), 
                     labels = c("Inferred trajectories" = "Inferred trajectories",
                                "True trajectories" = "True trajectories",
                                "train_obs" = "Observations for fitting", 
                                "test_obs" = "Observations for prediction"
                     ),
                     name = NULL) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape = FALSE,
         col =guide_legend(override.aes = list(linetype = c(0,0,1,1),
                                               shape = c(16,1,NA,NA),
                                               color = c("black" , "black",
                                                         line_color,
                                                         "red"))))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab(NULL) + ylab("R") + geom_vline(xintercept = train_day, linetype = "dashed")


p_4 <- ggplot() +
  geom_ribbon(data = df_D, aes(ymin = q.025, ymax = q.975, x = time, fill = is_train_D), alpha = 0.4) +
  geom_line(data = data_df_truth, aes(x = time, y = D, color = "True trajectories"), size = 3) +
  geom_line(data = df_D, aes(x = time, y = est, color = "Inferred trajectories"), size = 1.75) +
  geom_point(data = data_df_obs, aes(x = time, y = D, shape = is_train_obs, color = is_train_obs), size = 3) +
  scale_fill_manual(values = c("train" = train_fill, "test" = predict_fill), 
                    labels = c("train" = "95% credible interval in fitting", 
                               "test" = "95% credible interval in prediction"), 
                    name = NULL) +
  scale_shape_manual(values = c("Inferred trajectories" = NA, 
                                "train_obs" = 16, "test_obs" = 1), 
                     labels = c("Inferred trajectories" = "Inferred trajectories",
                                "train_obs" = "Observations for fitting", 
                                "test_obs" = "Observations for prediction"),
                     name = NULL) +
  scale_color_manual(values = c("Inferred trajectories" = line_color,
                                "True trajectories" = "red",
                                "train_obs" = "black", "test_obs" = "black"), 
                     labels = c("Inferred trajectories" = "Inferred trajectories",
                                "True trajectories" = "True trajectories",
                                "train_obs" = "Observations for fitting", 
                                "test_obs" = "Observations for prediction"
                     ),
                     name = NULL) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape = FALSE,
         col =guide_legend(override.aes = list(linetype = c(0,0,1,1),
                                               shape = c(16,1,NA,NA),
                                               color = c("black" , "black",
                                                         line_color,
                                                         "red"))))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab(NULL) + ylab("D") + geom_vline(xintercept = train_day, linetype = "dashed")


# Combine plots
p <- p_1 + p_2 + p_3 + p_4 + plot_layout(guides = "collect", ncol = 2) & 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.box = "horizontal") 

p
