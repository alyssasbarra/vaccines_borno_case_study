input_data <- fread("FILEPATH")
input_data[, age_cohort := original_year - year + 1]

nid_age_cov <- input_data[, .(wt_cov = weighted.mean(dpt3_cov/N, w = N * weight),
                              N = sum(N * weight)), 
							  by = c("svy_id", "age_cohort")]

baseline_cov <- subset(nid_age_cov, age_cohort == 1, select = c("wt_cov", "svy_id"))
setnames(baseline_cov, "wt_cov", "baseline_wt_cov")

nid_age_cov <- merge(nid_age_cov, baseline_cov, by = "svy_id")
nid_age_cov[, diff := wt_cov - baseline_wt_cov]

cohort_age_cov <- nid_age_cov[, .(diff = weighted.mean(diff, w = N)),
				   				by = "age_cohort"]

nid_age_cov$age_cohort <- factor(nid_age_cov$age_cohort)

gg <- ggplot() +
	  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
	  geom_line(data = nid_age_cov, 
             aes(x = age_cohort, y = diff, group = svy_id),
             alpha = 0.08) +
	  geom_line(data = cohort_age_cov, 
	            aes(x = age_cohort, y = diff), color = "red", size = 1, alpha = 0.5) +
	  theme_classic() +
	  labs(x = "Age Cohort", 
	       y = "DPT3 coverage difference\n(coverage in age cohort of interest - coverage in 12-23 month olds)") +
	  scale_x_discrete(breaks=c(1,2,3,4),
        labels=c("12-23 months", "24-35 months", "36-47 months", "48-59 months"),
        expand = c(0,0)) +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1))


png(file = "FILEPATH",
    width = 10,
    height = 6,
    units = "in",
    res = 300)

print(gg)
dev.off()