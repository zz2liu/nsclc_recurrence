# Plot Surv(time = Months, event = Status) ~ {strat} to PDF and show.
# Also save medianSurv and survRate each year as html files.
km_plot = function(data, strat, pval=T, 
                   out_prefix='out/', xlim=c(0, 120), ...) {
  f = as.formula(glue('Surv(time = Months, event = Status) ~ {strat}'))
  fit = surv_fit(f, data=data)
  
  smed = tbl_survfit(fit,
                     probs = 0.5,
                     label_header = "**Median survival (95% CI)**")
  smed %>% as_gt() %>% gtsave(glue("{out_prefix}_medianSurv.html"))
  
  srate = tbl_survfit(fit,
                      times = seq(12, 60, by=12),
                      label_header = '**Year {time/12}**')
  srate %>% as_gt() %>% gtsave(glue("{out_prefix}_survRate.html")) #.rtf 
  
  survp = ggsurvplot(fit #, data=data #but lose the risktable
                     , break.x.by=12
                     , pval=pval
                     , risk.table='abs_pct', fontsize=3
                     , surv.median.line='hv'
                     , xlim=xlim)
  pdf(glue("{out_prefix}_survPlot.pdf"), ...)
  print(survp, newpage=F)
  dev.off()
  survp
}

ggsurvfit_plot = function(data, strat, pval=T, 
                   out_prefix='out/', ...) {
  f = as.formula(glue('Surv(time = Months, event = Status) ~ {strat}'))
  fit = survfit2(f, data=data)
  
  smed = tbl_survfit(fit,
                     probs = 0.5,
                     label_header = "**Median survival (95% CI)**")
  smed %>% as_gt() %>% gtsave(glue("{out_prefix}_medianSurv.html"))
  
  srate = tbl_survfit(fit,
                      times = seq(12, 60, by=12),
                      label_header = '**Year {time/12}**')
  srate %>% as_gt() %>% gtsave(glue("{out_prefix}_survRate.html")) #.rtf 
  
  fit |>
    ggsurvfit(linewidth = 1) +
    add_censor_mark() +
    add_confidence_interval() +
    add_pvalue(size = 5) +
    add_risktable() +
    add_quantile(y_value = 0.6, color = "gray50", linewidth = 0.75) +
    scale_ggsurvfit(x_scales = list(breaks = seq(0, 144, by = 12))) +
    coord_cartesian(xlim = c(0, 120)) +
    labs(x='Months')
  ggsave(glue("{out_prefix}_survPlot.pdf"), ...)
}

