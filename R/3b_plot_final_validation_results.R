
# 
# selected_models |>
#   left_join(models_rs |> select(model_type, model_id, final_metrics))
# 
# prev_best_models |>
#   unnest(reg_cv_res)
# fit_models |>
#   filter(model_type == 'multinomial regression') |>
#   unnest(params) |>
#   ggplot(aes(penalty, mixture)) +
#   geom_text(aes(label = model_id)) +
#   scale_x_log10()
# 
# model_rs |>
#   select(model_type, model_id, final_metrics) |>
#   unnest(final_metrics) |>
#   filter(.metric == 'mcc') |>
#   ggplot(aes(x = .estimate, y = fct_rev(factor(model_id)), fill = model_type)) +
#   facet_wrap(~model_type) +
#   geom_col()

# 
# fitted_models |> 
#   filter(model_type == 'multinom_reg_glmnet' & model_id == 7) |> 
#   select(preds) |> unnest(cols = c(preds)) |> 
#   bind_cols(truth = test_prep$subfamily) |> 
#   my_metrics(truth = truth, estimate = .pred_class, na_rm = T)
#   
# # wtf, now model 7 is terrible? This is where I found issue arising from max_entropy grid creating models randomly each time -> ids did not match up when `models` created multiple times
# 
# 
