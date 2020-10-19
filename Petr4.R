# QUANDO A PETR4 VOLTA ?

setwd("C:/Users/lucmo/Downloads/Projetos_Programacao/R/Fincance/Petr4/")

graphics.off()
library(BatchGetSymbols)
library(tidyverse)
library(quantmod)
library(rugarch)
first_date <- '2016-01-01' 
last_date <- '2020-10-16' 
my_ticker <- 'PETR4.SA' 
series_name <- 'PETR4.SA' 

l_out <- BatchGetSymbols(tickers = my_ticker, 
                         first.date = first_date, 
                         last.date = last_date)

df_prices <- l_out$df.tickers %>%
  select(ref.date, ticker, price.adjusted) %>%
  mutate(log_ret = log(price.adjusted/dplyr::lag(price.adjusted) ),
         arim_ret = price.adjusted/dplyr::lag(price.adjusted) - 1,
         series_name = series_name) %>%
  na.omit() # remove all NA values

View(df_prices)

petr4.xts = xts(df_prices$log_ret,order.by =df_prices$ref.date )

chartSeries(petr4.xts, subset = "2019-10-18/2020-10-15")


#second part


do_arch_test(x = df_prices$log_ret, max_lag = 10) # teste LM

ar_lag <- 0 # lag used for ar term in mean equation (0 in paper)
ma_lag <- 0 # lag used for ma term in mean equation (0 in paper)
arch_lag <- 1 # lag in arch effect (1 in paper)
garch_lag <- 1 # lag in garch effect (1 in paper)
models_to_estimate <- c('sGARCH', 'eGARCH', 'gjrGARCH') # see rugarch manual for more
distribution_to_estimate <- 'norm' # distribution used in all models
df_grid <- expand_grid(ar_lag,
                       ma_lag,
                       arch_lag,
                       garch_lag,
                       models_to_estimate,
                       distribution_to_estimate)










estimate_garch <- function(ar_lag,
                           ma_lag,
                           arch_lag,
                           garch_lag,
                           models_to_estimate,
                           distribution_to_estimate) {
  
  message('Estimating ARMA(',ar_lag,',', ma_lag, ')', '-',
          models_to_estimate, '(', arch_lag, ',', garch_lag, ') ', 
          'dist = ', distribution_to_estimate)
  
  # estimate model
  my_spec <- ugarchspec(variance.model = list(model = models_to_estimate,
                                              garchOrder = c(arch_lag, 
                                                             garch_lag)),
                        mean.model = list(armaOrder = c(ar_lag,
                                                        ma_lag)), 
                        distribution.model = distribution_to_estimate)
  
  my_garch <- ugarchfit(spec = my_spec, data = df_prices$log_ret)
  
  return(my_garch)
}


l_args <- as.list(df_grid)
l_models <- pmap(.l = l_args, .f = estimate_garch)
l_models <- map(l_models, extract.rugarch, include.rsquared = FALSE)


l_models

#5


max_lag_AR <- 1 # used 1 in paper
max_lag_MA <- 1 # used 1 in paper
max_lag_ARCH <- 2 # used 2 in paper
max_lag_GARCH <- 1 # used 1 in paper
dist_to_use <- c('norm', 'std') # see rugarch::ugarchspec help for more
models_to_estimate <- c('sGARCH', 'eGARCH', 'gjrGARCH') # see rugarch::rugarchspec help for more



out <- find_best_arch_model(x = df_prices$log_ret, 
                            type_models = models_to_estimate,
                            dist_to_use = dist_to_use,
                            max_lag_AR = max_lag_AR,
                            max_lag_MA = max_lag_MA,
                            max_lag_ARCH = max_lag_ARCH,
                            max_lag_GARCH = max_lag_GARCH)

# get table with estimation results
tab_out <- out$tab_out

# pivot table to long format (better for plotting)
df_long <- tidyr::pivot_longer(data = tab_out %>%
                                 select(model_name,
                                        type_model,
                                        type_dist,
                                        AIC, BIC),  cols = c('AIC', 'BIC'))

models_names <- unique(df_long$model_name)
best_models <- c(tab_out$model_name[which.min(tab_out$AIC)],
                 tab_out$model_name[which.min(tab_out$BIC)])
df_long <- df_long %>%
  mutate(order_model = if_else(model_name %in% best_models, 'Best Model', 'Not Best Model') ) %>%
  na.omit()

# make table with best models
df_best_models <- df_long %>%
  group_by(name) %>%
  summarise(model_name = model_name[which.min(value)],
            value = value[which.min(value)],
            type_model = type_model[which.min(value)])

# plot results
ggplot(df_long %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  theme_bw(base_family = "TT Times New Roman") + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting Garch Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")



# estimate best garch model by BIC (used in next section)
best_spec = ugarchspec(variance.model = list(model =  out$best_bic$type_model, 
                                             garchOrder = c(out$best_bic$lag_arch,
                                                            out$best_bic$lag_garch)),
                       mean.model = list(armaOrder = c(out$best_bic$lag_ar, 
                                                       out$best_bic$lag_ma)),
                       distribution = 'std')

my_best_garch <- ugarchfit(spec = best_spec, 
                           data = df_prices$log_ret)



set.seed(1982) # fix seed for simulations (20200315 replicates the paper's results)
n_sim <- 5000


n_days_ahead <- 2*365 


graphics.off()
my_garch<-my_best_garch

df_sim <- do_sim(n_sim = n_sim, 
                n_t = n_days_ahead, 
                my_garch, 
                df_prices = df_prices)


df_sim %>% View()
df_sim %>% glimpse()


# calculate probabilities of reaching peak value
tab_prob <- df_sim %>%
  group_by(ref_date) %>%
  summarise(prob = mean(sim_price > max(df_prices$price.adjusted)))

n_years_back <- 4
df_prices_temp <- df_prices %>%
  dplyr::filter(ref.date > max(ref.date) - n_years_back*365)
my_garch_name <- toupper(as.character(my_garch@model$modeldesc$vmodel))


ggplot() + 
  geom_line(data = df_prices_temp, 
            aes(x = ref.date, y = price.adjusted), color = 'black', size = 0.75)  + 
  geom_line(data = df_sim, 
            aes(x = ref_date, 
                y = sim_price, 
                group = i_sim),
            color = 'grey', 
            size = 0.05,
            alpha = 0.015) + 
  theme_bw(base_family = "TT Times New Roman") + 
  geom_hline(yintercept = max(df_prices_temp$price.adjusted),color='red') + 
  labs(title = paste0('Projeções de Preço ', series_name),
       subtitle = paste0('Total de ', n_sim, ' simulações de preço baseadas no modelo ',
                         my_garch_name, 
                         'selecionado por BIC'),
       caption = 'Dados do Yahoo Finance',
       x = '',
       y = 'Preço') + 
  ylim(c(0.75*min(df_prices_temp$price.adjusted), 
         1.25*max(df_prices_temp$price.adjusted))) + 
  xlim(c(max(df_prices_temp$ref.date) - n_years_back*365,
         max(df_prices_temp$ref.date) + 2*365) )



