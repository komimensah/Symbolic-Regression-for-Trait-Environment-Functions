rm(list=ls()); gc()

library(tidyverse)
library(reticulate)
library(minpack.lm)
setwd('/Users/kagboka/Desktop/Stomoxys calcitrans/Ecological_modelling/Figures/')
# ============================================================
# USE PYTHON ENVIRONMENT
# ============================================================

use_virtualenv("r-pysr-env", required = TRUE)
pysr <- import("pysr")

# ============================================================
# INPUT DATA
# ============================================================

adult <- tribble(
  ~Diet, ~Temp, ~Days,
  "100%",15,60,
  "100%",20,25,
  "100%",25,12,
  "100%",30,8,
  "100%",35,10,
  "50%",15,70,
  "50%",20,30,
  "50%",25,15,
  "50%",30,10,
  "50%",35,12,
  "25%",15,80,
  "25%",20,35,
  "25%",25,18,
  "25%",30,12,
  "25%",35,15,
  "12%",15,90,
  "12%",20,40,
  "12%",25,20,
  "12%",30,14,
  "12%",35,16
) %>% mutate(Stage="Adult")

pupal <- tribble(
  ~Diet, ~Temp, ~Days,
  "100%",15,60,
  "100%",20,25,
  "100%",25,12,
  "100%",30,8,
  "100%",35,10,
  "50%",15,65,
  "50%",20,30,
  "50%",25,15,
  "50%",30,10,
  "50%",35,11,
  "25%",15,85,
  "25%",20,35,
  "25%",25,18,
  "25%",30,12,
  "25%",35,12,
  "12%",15,95,
  "12%",20,40,
  "12%",25,20,
  "12%",30,14,
  "12%",35,13
) %>% mutate(Stage="Pupal")

df <- bind_rows(adult, pupal) %>%
  mutate(rate = 1/Days)

# ============================================================
# COMMON METRIC FUNCTION
# ============================================================

evaluate_model <- function(y_true, y_pred, k){
  residuals <- y_true - y_pred
  rss <- sum(residuals^2)
  tss <- sum((y_true - mean(y_true))^2)
  r2 <- 1 - rss/tss
  rmse <- sqrt(mean(residuals^2))
  n <- length(y_true)
  aic <- n*log(rss/n) + 2*k
  return(list(R2=r2, RMSE=rmse, AIC=aic))
}

# --------------------------------------------------
# SYMBOLIC REGRESSION (CONSTRAINED)
# --------------------------------------------------
fit_symbolic_full <- function(dat){
  
  X <- matrix(dat$Temp, ncol = 1)
  y <- dat$rate
  
  # --- Biological constraint #1: positivity ---
  # Fit on log-rate to enforce positive predictions after exp()
  eps <- 1e-8
  y_log <- log(pmax(y, eps))
  
  model <- pysr$PySRRegressor(
    niterations = as.integer(400),
    npop        = as.integer(60),
    
    # --- Biological constraint #2: remove division (prevents poles/spikes) ---
    binary_operators = list("+", "-", "*"),
    
    # --- Biological constraint #3: safe unary operators ---
    # sqrt is ok only if used safely; exp helps represent accelerating rates
    unary_operators  = list("sqrt", "exp"),
    
    # Encourage simpler equations
    # (keeps comparison fair + improves interpretability)
    parsimony = 1e-2,
    
    verbosity = 0
  )
  
  model$fit(X, y_log)
  
  # Predictions on training temps
  yhat_log <- as.numeric(model$predict(X))
  yhat <- exp(yhat_log)   # back-transform => always positive
  
  metrics <- evaluate_model(y, yhat, k = 3)
  
  list(
    predictions = yhat,
    equation_logscale = model$get_best()$equation,  # equation for log(rate)
    metrics = metrics,
    model = model
  )
}

# --------------------------------------------------
# APPLY TO ALL DIET × STAGE
# --------------------------------------------------
symbolic_results <- list()
symbolic_plot_df <- data.frame()

temp_grid <- seq(15, 35, length.out = 200)

for(st in unique(df$Stage)){
  for(d in unique(df$Diet)){
    
    subset_dat <- df %>% filter(Stage == st, Diet == d)
    
    fit <- fit_symbolic_full(subset_dat)
    
    # store table row
    symbolic_results[[paste(st, d, sep = "_")]] <- tibble(
      Stage = st,
      Diet  = d,
      Model = "Symbolic (constrained)",
      Equation = paste0("log(rate) = ", fit$equation_logscale, " ; rate = exp(.)"),
      R2   = fit$metrics$R2,
      RMSE = fit$metrics$RMSE,
      AIC  = fit$metrics$AIC
    )
    
    # smooth curve on grid (still positive)
    Xnew <- matrix(temp_grid, ncol = 1)
    preds_log <- as.numeric(fit$model$predict(Xnew))
    preds <- exp(preds_log)
    
    symbolic_plot_df <- rbind(
      symbolic_plot_df,
      data.frame(Stage = st, Diet = d, Temp = temp_grid, rate = preds)
    )
  }
}

symbolic_table <- bind_rows(symbolic_results)

print(symbolic_table)

# --------------------------------------------------
# SAVE SYMBOLIC FIGURE (2 PANELS: Adult | Pupal)
# --------------------------------------------------
symbolic_plot <- ggplot() +
  geom_point(data = df, aes(Temp, rate, color = Diet), size = 3) +
  geom_line(data = symbolic_plot_df, aes(Temp, rate, color = Diet), linewidth = 1.2) +
  facet_wrap(~Stage, nrow = 1) +
  theme_classic(base_size = 14) +
  labs(
    title = "Symbolic Regression (biologically constrained: positive, no singularities)",
    x = "Temperature (°C)",
    y = "Development rate (1/days)",
    color = "Diet"
  )

ggsave("Symbolic_models_constrained.png", symbolic_plot, width = 10, height = 5, dpi = 300)


# ============================================================
# CLASSICAL CANONICAL DEVELOPMENT MODELS
# ============================================================
library(minpack.lm)
library(dplyr)
library(ggplot2)

# --------------------------------------------------
# MODEL EVALUATION FUNCTION
# --------------------------------------------------
evaluate_model <- function(y_true, y_pred, k){
  residuals <- y_true - y_pred
  rss <- sum(residuals^2)
  tss <- sum((y_true - mean(y_true))^2)
  r2 <- 1 - rss/tss
  rmse <- sqrt(mean(residuals^2))
  n <- length(y_true)
  aic <- n*log(rss/n) + 2*k
  list(R2=r2, RMSE=rmse, AIC=aic)
}

# --------------------------------------------------
# FIT ALL CANONICAL MODELS
# --------------------------------------------------
fit_classical <- function(dat){
  
  results <- list()
  
  try_model <- function(name, formula, start, k){
    
    fit <- tryCatch(
      nlsLM(formula,
            data = dat,
            start = start,
            control = nls.lm.control(maxiter=1000)),
      error = function(e) NULL
    )
    
    if(!is.null(fit)){
      pred <- predict(fit)
      metrics <- evaluate_model(dat$rate, pred, k)
      
      results[[name]] <<- list(
        params = coef(fit),
        equation = deparse(formula),
        metrics = metrics
      )
    }
  }
  
  # Quadratic
  try_model("Quadratic",
            rate ~ a + b*Temp + c*Temp^2,
            list(a=0.01,b=0.001,c=-0.00001),3)
  
  # Briere
  try_model("Briere",
            rate ~ a*Temp*(Temp-Tmin)*sqrt(Tmax-Temp),
            list(a=0.0001,Tmin=10,Tmax=40),3)
  
  # Logan
  try_model("Logan",
            rate ~ psi*(exp(rho*Temp) - exp(rho*Tmax - (Tmax-Temp)/Delta)),
            list(psi=0.1,rho=0.01,Tmax=40,Delta=1),4)
  
  # Lactin
  try_model("Lactin",
            rate ~ exp(rho*Temp) - exp(rho*Tmax - (Tmax-Temp)/Delta) + lambda,
            list(rho=0.01,Tmax=40,Delta=1,lambda=0.01),4)
  
  # Ratkowsky
  try_model("Ratkowsky",
            rate ~ (b*(Temp-Tmin))^2,
            list(b=0.01,Tmin=10),2)

  
  results
}

# --------------------------------------------------
# FIT PER Stage × Diet AND STORE RESULTS
# --------------------------------------------------

temp_grid <- seq(15,35,length.out=200)

classical_results_list <- list()
classical_plot_df_clean <- data.frame()

for(st in unique(df$Stage)){
  for(d in unique(df$Diet)){
    
    subset_dat <- df %>% filter(Stage==st, Diet==d)
    fits <- fit_classical(subset_dat)
    
    if(length(fits)==0) next
    
    # Build summary table for this group
    group_table <- bind_rows(
      lapply(names(fits), function(m){
        data.frame(
          Stage=st,
          Diet=d,
          Model=m,
          Equation=fits[[m]]$equation,
          R2=fits[[m]]$metrics$R2,
          RMSE=fits[[m]]$metrics$RMSE,
          AIC=fits[[m]]$metrics$AIC
        )
      })
    )
    
    # Select best model by AIC
    best_row <- group_table %>% slice_min(AIC,n=1,with_ties=FALSE)
    
    classical_results_list[[paste(st,d,sep="_")]] <- best_row
    
    # Reconstruct predictions for best model
    best_model <- best_row$Model
    params <- fits[[best_model]]$params
    
    ypred <- switch(best_model,
                    
                    "Quadratic" =
                      params["a"] + params["b"]*temp_grid + params["c"]*temp_grid^2,
                    
                    "Briere" =
                      params["a"]*temp_grid*
                      (temp_grid-params["Tmin"])*
                      sqrt(params["Tmax"]-temp_grid),
                    
                    "Logan" =
                      params["psi"]*
                      (exp(params["rho"]*temp_grid) -
                         exp(params["rho"]*params["Tmax"] -
                               (params["Tmax"]-temp_grid)/params["Delta"])),
                    
                    "Lactin" =
                      exp(params["rho"]*temp_grid) -
                      exp(params["rho"]*params["Tmax"] -
                            (params["Tmax"]-temp_grid)/params["Delta"]) +
                      params["lambda"],
                    
                    "Ratkowsky" =
                      (params["b"]*(temp_grid-params["Tmin"]))^2
    )
    
    classical_plot_df_clean <- rbind(
      classical_plot_df_clean,
      data.frame(Stage=st,Diet=d,Temp=temp_grid,rate=ypred)
    )
  }
}

# Final comparison table
classical_table <- bind_rows(classical_results_list)

print(classical_table)

# --------------------------------------------------
# CLEAN CLASSICAL FIGURE (MATCHES SYMBOLIC STYLE)
# --------------------------------------------------

classical_plot_clean <- ggplot() +
  geom_point(data=df,
             aes(Temp,rate,color=Diet),
             size=3) +
  geom_line(data=classical_plot_df_clean,
            aes(Temp,rate,color=Diet),
            linewidth=1.2) +
  facet_wrap(~Stage,nrow=1) +
  theme_classic(base_size=14) +
  labs(
    title="Classical Development Models (Best per Diet × Stage)",
    x="Temperature (°C)",
    y="Development rate (1/days)",
    color="Diet"
  )

ggsave("Classical_models_clean.png",
       classical_plot_clean,
       width=10,height=5,dpi=300)
# ============================================================
# CONSOLIDATED COMPARISON TABLE
# ============================================================

symbolic_table$ModelType <- "Symbolic"
classical_table$ModelType <- "Classical"

comparison_table <- bind_rows(
  symbolic_table %>% mutate(Model="Symbolic") %>%
    select(Stage,Diet,ModelType,Model,Equation,R2,RMSE,AIC),
  classical_table %>%
    select(Stage,Diet,ModelType,Model,Equation,R2,RMSE,AIC)
)

write.csv(comparison_table,
          "Model_Comparison_Table.csv",
          row.names=FALSE)

print(comparison_table)
