
rm(list=ls()); gc()

library(tidyverse)
library(reticulate)
library(minpack.lm)

setwd('/Users/kagboka/Desktop/Stomoxys calcitrans/Ecological_modelling/Figures/')

# ============================================================
# PYTHON ENV
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
  list(R2=r2, RMSE=rmse, AIC=aic)
}

# ============================================================
# SYMBOLIC REGRESSION
# ============================================================

fit_symbolic_full <- function(dat){
  
  X <- matrix(dat$Temp, ncol = 1)
  y <- dat$rate
  
  eps <- 1e-8
  y_log <- log(pmax(y, eps))
  
  model <- pysr$PySRRegressor(
    niterations = as.integer(400),
    npop = as.integer(60),
    binary_operators = list("+","-","*"),
    unary_operators = list("sqrt","exp"),
    parsimony = 1e-2,
    verbosity = 0
  )
  
  model$fit(X, y_log)
  
  yhat_log <- as.numeric(model$predict(X))
  yhat <- exp(yhat_log)
  
  metrics <- evaluate_model(y, yhat, k=3)
  
  list(
    model = model,
    equation_log = model$get_best()$equation,
    metrics = metrics
  )
}

symbolic_results <- list()
symbolic_model_store <- list()

for(st in unique(df$Stage)){
  
  symbolic_model_store[[st]] <- list()
  
  for(d in unique(df$Diet)){
    
    subdat <- df %>% filter(Stage==st, Diet==d)
    fit <- fit_symbolic_full(subdat)
    
    symbolic_model_store[[st]][[d]] <- fit$model
    
    symbolic_results[[paste(st,d,sep="_")]] <- tibble(
      Stage=st,
      Diet=d,
      Model="Symbolic",
      Equation=paste0("log(rate)=",fit$equation_log),
      R2=fit$metrics$R2,
      RMSE=fit$metrics$RMSE,
      AIC=fit$metrics$AIC
    )
  }
}

symbolic_table <- bind_rows(symbolic_results)
print(symbolic_table)

# ============================================================
# CLASSICAL MODELS
# ============================================================

fit_classical <- function(dat){
  
  results <- list()
  
  try_model <- function(name, formula, start, k){
    
    fit <- tryCatch(
      nlsLM(formula,
            data=dat,
            start=start,
            control=nls.lm.control(maxiter=1000)),
      error=function(e) NULL
    )
    
    if(!is.null(fit)){
      pred <- predict(fit)
      metrics <- evaluate_model(dat$rate, pred, k)
      results[[name]] <<- list(params=coef(fit),
                               equation=deparse(formula),
                               metrics=metrics)
    }
  }
  
  try_model("Quadratic",
            rate ~ a + b*Temp + c*Temp^2,
            list(a=0.01,b=0.001,c=-0.00001),3)
  
  try_model("Briere",
            rate ~ a*Temp*(Temp-Tmin)*sqrt(Tmax-Temp),
            list(a=0.0001,Tmin=10,Tmax=40),3)
  
  try_model("Logan",
            rate ~ psi*(exp(rho*Temp) - exp(rho*Tmax - (Tmax-Temp)/Delta)),
            list(psi=0.1,rho=0.01,Tmax=40,Delta=1),4)
  
  try_model("Lactin",
            rate ~ exp(rho*Temp) - exp(rho*Tmax - (Tmax-Temp)/Delta) + lambda,
            list(rho=0.01,Tmax=40,Delta=1,lambda=0.01),4)
  
  try_model("Ratkowsky",
            rate ~ (b*(Temp-Tmin))^2,
            list(b=0.01,Tmin=10),2)
  
  results
}

classical_best_store <- list()
classical_results <- list()

for(st in unique(df$Stage)){
  
  classical_best_store[[st]] <- list()
  
  for(d in unique(df$Diet)){
    
    subdat <- df %>% filter(Stage==st, Diet==d)
    fits <- fit_classical(subdat)
    if(length(fits)==0) next
    
    table_tmp <- bind_rows(
      lapply(names(fits), function(m){
        data.frame(
          Model=m,
          R2=fits[[m]]$metrics$R2,
          RMSE=fits[[m]]$metrics$RMSE,
          AIC=fits[[m]]$metrics$AIC
        )
      })
    )
    
    best <- table_tmp %>% slice_min(AIC,n=1,with_ties=FALSE)
    best_model <- best$Model
    
    classical_best_store[[st]][[d]] <- list(
      model=best_model,
      params=fits[[best_model]]$params
    )
    
    classical_results[[paste(st,d,sep="_")]] <- tibble(
      Stage=st,
      Diet=d,
      Model=best_model,
      Equation=fits[[best_model]]$equation,
      R2=fits[[best_model]]$metrics$R2,
      RMSE=fits[[best_model]]$metrics$RMSE,
      AIC=fits[[best_model]]$metrics$AIC
    )
  }
}

classical_table <- bind_rows(classical_results)
print(classical_table)

# ============================================================
# CONSOLIDATED TABLE
# ============================================================

symbolic_table$ModelType <- "Symbolic"
classical_table$ModelType <- "Classical"

comparison_table <- bind_rows(
  symbolic_table,
  classical_table
)
print(comparison_table)

library(terra)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(dplyr)

# --------------------------------------------------
# 1) LOAD YOUR REAL TEMPERATURE RASTER
# --------------------------------------------------

T_ne <- rast('/Users/kagboka/Desktop/Lesh2.0/Temperature.tif')

# Nebraska boundary
usa <- ne_states(country="United States of America", returnclass="sf")
nebraska <- usa %>% filter(name=="Nebraska")
nebraska_vect <- vect(nebraska)

# Crop + mask to Nebraska
T_ne <- crop(T_ne, nebraska_vect)
T_ne <- mask(T_ne, nebraska_vect)

temp_values <- values(T_ne)

# Valid cells only
valid_idx <- which(!is.na(temp_values))

# Optional: clamp to training range (highly recommended)
temp_values[temp_values < 15] <- 15
temp_values[temp_values > 35] <- 35

# --------------------------------------------------
# 2) PREDICT SYMBOLIC & CLASSICAL
# --------------------------------------------------

symbolic_rasters  <- list()
classical_rasters <- list()

for(st in names(symbolic_model_store)){
  
  symbolic_rasters[[st]]  <- list()
  classical_rasters[[st]] <- list()
  
  for(d in names(symbolic_model_store[[st]])){
    
    # ================= SYMBOLIC =================
    
    sym_model <- symbolic_model_store[[st]][[d]]
    sym_pred <- rep(NA_real_, length(temp_values))
    
    if(length(valid_idx) > 0){
      
      Xnew <- matrix(temp_values[valid_idx], ncol=1)
      sym_pred_log <- as.numeric(sym_model$predict(Xnew))
      sym_vals <- exp(sym_pred_log)
      
      # Remove invalid numerical results
      sym_vals[!is.finite(sym_vals)] <- NA
      
      sym_pred[valid_idx] <- sym_vals
    }
    
    sym_r <- T_ne
    values(sym_r) <- sym_pred
    symbolic_rasters[[st]][[d]] <- sym_r
    
    
    # ================= CLASSICAL =================
    
    cls_info   <- classical_best_store[[st]][[d]]
    best_model <- cls_info$model
    params     <- cls_info$params
    
    cls_pred <- rep(NA_real_, length(temp_values))
    
    if(length(valid_idx) > 0){
      
      Tv <- temp_values[valid_idx]
      
      pred_vals <- switch(best_model,
                          
                          "Quadratic" =
                            params["a"] + params["b"]*Tv + params["c"]*Tv^2,
                          
                          "Briere" =
                            params["a"]*Tv*(Tv-params["Tmin"])*
                            sqrt(pmax(params["Tmax"]-Tv,0)),
                          
                          "Logan" =
                            params["psi"]*
                            (exp(params["rho"]*Tv) -
                               exp(params["rho"]*params["Tmax"] -
                                     (params["Tmax"]-Tv)/params["Delta"])),
                          
                          "Lactin" =
                            exp(params["rho"]*Tv) -
                            exp(params["rho"]*params["Tmax"] -
                                  (params["Tmax"]-Tv)/params["Delta"]) +
                            params["lambda"],
                          
                          "Ratkowsky" =
                            (params["b"]*(Tv-params["Tmin"]))^2
      )
      
      # Remove invalid numerical results
      pred_vals[!is.finite(pred_vals)] <- NA
      
      cls_pred[valid_idx] <- pred_vals
    }
    
    cls_r <- T_ne
    values(cls_r) <- cls_pred
    classical_rasters[[st]][[d]] <- cls_r
  }
}

# --------------------------------------------------
# 3) BUILD DELTA RASTERS
# --------------------------------------------------

delta_rasters <- list()

for(st in names(symbolic_rasters)){
  
  delta_rasters[[st]] <- list()
  
  for(d in names(symbolic_rasters[[st]])){
    
    sym_r <- symbolic_rasters[[st]][[d]]
    cls_r <- classical_rasters[[st]][[d]]
    
    if(!is.null(sym_r) && !is.null(cls_r)){
      delta_rasters[[st]][[d]] <- sym_r - cls_r
    }
  }
}

# --------------------------------------------------
# 4) CONVERT TO DATAFRAME (SAFE)
# --------------------------------------------------

map_df <- data.frame()

for(st in names(delta_rasters)){
  for(d in names(delta_rasters[[st]])){
    
    r_obj <- delta_rasters[[st]][[d]]
    if(is.null(r_obj)) next
    
    tmp <- as.data.frame(r_obj, xy=TRUE)
    tmp <- tmp[!is.na(tmp[,3]), ]
    
    if(nrow(tmp) == 0) next
    
    names(tmp)[3] <- "value"
    tmp$Stage <- st
    tmp$Diet  <- d
    
    map_df <- rbind(map_df, tmp)
  }
}

map_df$Stage <- factor(map_df$Stage, levels=c("Adult","Pupal"))
map_df$Diet  <- factor(map_df$Diet, levels=c("12%","25%","50%","100%"))

# --------------------------------------------------
# 5) PLOT DELTA MAP
# --------------------------------------------------

delta_plot <- ggplot(map_df, aes(x, y, fill = value)) +
  geom_raster() +
  geom_sf(
    data = nebraska,
    fill = NA,
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  facet_grid(Stage ~ Diet) +
  coord_sf(expand = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4)
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    fill = "Δ rate",
    title = "Nebraska: Symbolic − Classical Development Rate"
  )

# Save at 300 dpi
ggsave(
  filename = "Nebraska_Delta_Symbolic_vs_Classical.png",
  plot = delta_plot,
  width = 10,
  height = 6,
  dpi = 300
)
# --------------------------------------------------
# LOAD SURROUNDING STATES ONLY
# --------------------------------------------------

usa_all <- ne_states(country = "United States of America", returnclass = "sf")

# Keep only Nebraska + neighbors
neighbors <- usa_all %>%
  filter(name %in% c("Nebraska",
                     "Kansas",
                     "Iowa",
                     "South Dakota",
                     "Wyoming",
                     "Colorado",
                     "Missouri"))

# --------------------------------------------------
# PLOT WITH LOCAL CONTEXT
# --------------------------------------------------

delta_plot <- ggplot() +
  
  # Neighboring states (light background)
  geom_sf(
    data = neighbors,
    fill = "grey95",
    color = "grey70",
    linewidth = 0.3
  ) +
  
  # Raster layer
  geom_raster(
    data = map_df,
    aes(x, y, fill = value)
  ) +
  
  # Emphasize Nebraska
  geom_sf(
    data = nebraska,
    fill = NA,
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  
  facet_grid(Stage ~ Diet) +
  coord_sf(expand = FALSE) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4)
  ) +
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  
  labs(
    fill = "Δ rate",
    title = "Nebraska (USA): Symbolic − Classical Development Rate"
  )

# Save
ggsave(
  "Nebraska_Delta_LocalContext.tiff",
  delta_plot,
  width = 11,
  height = 6,
  dpi = 300,
  compression = "lzw"
)
