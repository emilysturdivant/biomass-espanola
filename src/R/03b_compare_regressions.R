# Compare model variants

library('tidyverse')

# Initialize
code <- 'HV_nu'
modeling_dir <- file.path('data/modeling', code)

# Read model results as list of tibbles
fps <- list.files(modeling_dir, 'model_results_.*\\.csv', 
                  recursive = TRUE, 
                  full.names = TRUE)
variants <- fps %>% str_extract(str_glue('(?<={code}/).*(?=/cal)'))
out <- fps %>% map(read_csv)

# Format
out2 <- out %>% 
  map_dfr(pivot_wider, names_from = c(name, type)) %>% 
  select(Intercept = estimate_int_train, 
                Slope = estimate_g0_train, 
                Intcpt_95CI = ci_int_train, 
                Slope_95CI = ci_slope_train,
                Adj_r2 = adj.r.squared_train, 
                starts_with('rmse'), 
                starts_with('mae'),
                starts_with('Bias')
                ) %>% 
  mutate(across(where(is.numeric), ~ signif(.x, digits = 4))) %>% 
  bind_cols(tibble(Variant = variants), .)

# Save
report_fp <- file.path('data/reports', str_c('model_results_by_variant_', code, '.csv'))
out2 %>% write_csv(report_fp)

# Look
var_sort <- c('cappt2_conserv13', 'cappt2_conserv13_mean5', 'med5', 'lee11s10', 'maskLU_lee11s10_LCinterp')
out3 <- out2 %>% 
  mutate(Variant = factor(Variant, levels = var_sort)) %>% 
  arrange(Variant) %>% 
  pivot_longer(cols = Adj_r2:BiasSD_test)

out2 %>% 
  select(Variant, contains('SD')) %>% 
  pivot_longer(cols = contains('SD'), values_to = 'SD') %>% 
  mutate(name = str_extract(name, '.*(?=SD)'))

out2 %>% 
  select(Variant, contains('_test')) %>%
  pivot_longer(cols = contains('_test'), values_to = 'value') %>% 
  mutate(name = str_extract(name, '.*(?=_test)'))

out3 %>% 
  filter(str_detect(name, 'SD'))
  select(contains('SD'))
  

out3 %>% ggplot(aes(y=value, x=Variant)) +
  geom_point() +
  facet_wrap(vars(name), scales = 'free')
