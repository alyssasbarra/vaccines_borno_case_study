# Source directly from within main script

indicator <- "dpt3_cov"

library(ggplot2)
library(data.table)
library(tidyr)

in_dir <- paste0("FILEPATH")

out_dir <- paste0("FILEPATH")

ad2_df <- fread(paste0("FILEPATH"))

#Subset to years of interest
ad2_df <- subset(ad2_df, year == max(year_list) | year == min(year_list))

# Convenience
ad2_df[year==max(year_list), year := 1]
ad2_df[year==min(year_list), year := 0]

# Convert to wide
ad2_df <- ad2_df %>% 
	gather(metric, value, -(ADM0_CODE:year)) %>% 
	unite(metric_year, metric, year) %>%
	spread(metric_year, value) %>%
	as.data.table

# Grab populations
load(paste0(in_dir, indicator, "_raked_admin_draws_eb_bin0_0.RData"))
pops <- subset(admin_2, year == max(year_list), select = c(ADM2_CODE, pop))
ad2_df <- merge(ad2_df, pops, by = "ADM2_CODE", all.x= T, all.y = F)

# Figure out which had significant increases or decreases
decreases <- fread(paste0(in_dir, "/number_plugging/admin_2_",
                          max(year_list), " - ", min(year_list), "_decline.csv"))

increases <- fread(paste0(in_dir, "/number_plugging/admin_2_",
                          max(year_list), " - ", min(year_list), "_increase.csv"))

sig_ad2s <- as.numeric(unique(c(increases$ADM2_CODE, decreases$ADM2_CODE)))

ad2_df[, sig := ifelse(ADM2_CODE %in% sig_ad2s, T, F)]

# Get SDI
source('FILEPATH')
gbd <- get_covariate_estimates(covariate_id = 881)
gbd <- subset(gbd, year_id %in% c(min(year_list), max(year_list)))

# Get into shape for merge
gaul_to_loc_id <- fread("FILEPATH")
setnames(gaul_to_loc_id, "GAUL_CODE", "ADM0_CODE")
gaul_to_loc_id <- subset(gaul_to_loc_id, select = c("ADM0_CODE", "loc_id"))
setnames(gaul_to_loc_id, "loc_id", "location_id")
gbd <- merge(gbd, gaul_to_loc_id, by = "location_id")
gbd <- subset(gbd, ADM0_CODE %in% unique(ad2_df$ADM0_CODE))
gbd <- subset(gbd, select = c("ADM0_CODE", "year_id", "mean_value"))
setnames(gbd, "mean_value", "sdi")

gbd[, year_id := paste0("sdi_", year_id)]
gbd <- gbd %>% spread(year_id, sdi) %>% as.data.table
gbd[, sdi_change := get(paste0("sdi_", max(year_list))) - get(paste0("sdi_", min(year_list)))]

ad2_df <- merge(ad2_df, gbd, all.x=T, by = "ADM0_CODE")

ad2_df[, sdi_cut := cut(get(paste0("sdi_", max(year_list))), breaks = seq(0,1,by=0.2))]

# Drop Ma'tan al Sarra
ad2_df<-subset(ad2_df, ADM0_NAME != "Ma'tan al-Sarra")

# Add in some probabilities of > 80% in 2015
ad2_80p_target <- fread(paste0("FILEPATH"))
ad2_80p_target <- subset(ad2_80p_target, year == max(year_list),
                         select = c("ADM2_CODE", "p_above"))

ad2_df <- merge(ad2_df, ad2_80p_target, by= "ADM2_CODE")
ad2_df[, p_above_cut := cut(p_above, breaks = seq(0,1,by=0.1))]

# For faceted plots
ad2_df[ADM0_NAME == "Democratic Republic of the Congo", ADM0_NAME := "DRC"]
ad2_df[ADM0_NAME == "Sao Tome and Principe", ADM0_NAME := "STP"]
ad2_df[ADM0_NAME == "United Republic of Tanzania", ADM0_NAME := "Tanzania"]
ad2_df[ADM0_NAME == "Central African Republic", ADM0_NAME := "CAR"]

# Use carto color scheme?
carto_colors <- c("#5F4690","#38A6A5","#73AF48","#E17C05","#94346E","#6F4070","#994E95")

carto_diverging <- c("#009392","#39b185","#9ccb86","#e9e29c","#eeb479","#e88471","#cf597e")

big_colors <- c("#d34a53","#36dee6","#862b1b","#45bc8d","#d74d82","#64b764","#ad74d6","#92b440",
				"#5f6ed3","#caa331","#573586","#c2a957","#628bd5","#c8772c","#cd76c7","#627123",
				"#a54381","#c6673c","#ac455f","#d17160")

png(filename=paste0(out_dir, "ad2_", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=300)

gg <- ggplot() +
				geom_point(data = ad2_df,
				           aes(x = mean_0, xmin = lower_0, xmax = upper_0,
				               y = mean_1, ymin = lower_1, ymax = upper_1,
				               size = pop, color = p_above, shape = sig),
				           alpha = 0.6) +
				geom_abline(alpha = 0.6) +
				geom_abline(slope = 0, intercept = 0.8, 
				            color = "black", alpha = 0.6,
				            linetype = "dotted") +
			#	geom_errorbar(alpha = 0.01) +
			#	geom_errorbarh(alpha = 0.0) +
				scale_size_area() +
			#	scale_alpha_discrete(range = c(0.2, 0.8)) +
			#	scale_color_manual(values = carto_colors) +
			#	scale_colour_distiller(palette = "Reds") +
				scale_color_gradientn(colors = rev(carto_diverging)) +
			#	facet_wrap(~ADM0_NAME) +
				theme_classic() +
				scale_shape_manual(values = c(1,19)) +
				labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
				     y = paste0("Mean DPT3 Coverage: ", max(year_list)),
				     size = "Number of children \n< 5 years of age",
				     color = paste0("Probability of \nDPT3 coverage > 80% \n(", max(year_list), ")"),
				     shape = "High (>95%) probability \nof true change 2000-2015") +
				scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
				scale_y_continuous(expand = c(0, 0), limits = c(0,1))

print(gg)

dev.off()

# Facet this one

png(filename=paste0(out_dir, "ad2_by_adm0_", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=12, 
    height=14, 
    pointsize=12, 
    res=300)

gg <- ggplot() +
				geom_point(data = ad2_df,
				           aes(x = mean_0, xmin = lower_0, xmax = upper_0,
				               y = mean_1, ymin = lower_1, ymax = upper_1,
				               size = pop, color = p_above, shape = sig),
				           alpha = 0.6) +
				geom_abline(alpha = 0.6) +
				geom_abline(slope = 0, intercept = 0.8, 
				            color = "black", alpha = 0.6,
				            linetype = "dotted") +
			#	geom_errorbar(alpha = 0.01) +
			#	geom_errorbarh(alpha = 0.0) +
				scale_size_area() +
			#	scale_alpha_discrete(range = c(0.2, 0.8)) +
			#	scale_color_manual(values = carto_colors) +
			#	scale_colour_distiller(palette = "Reds") +
				scale_color_gradientn(colors = rev(carto_diverging), limits = c(0,1)) +
				facet_wrap(~ADM0_NAME, ncol = 6) +
				theme_minimal() +
				theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
				scale_shape_manual(values = c(1,19)) +
				labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
				     y = paste0("Mean DPT3 Coverage: ", max(year_list)),
				     size = "Number of children \n< 5 years of age",
				     color = paste0("Probability of \nDPT3 coverage > 80% \n(", max(year_list), ")"),
				     shape = "High (>95%) probability \nof true change 2000-2015") +
				scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
				scale_y_continuous(expand = c(0, 0), limits = c(0,1))

print(gg)

dev.off()

country_list <- unique(ad2_df$ADM0_NAME)
country_list <- list(country_list[1:18],
                     country_list[19:36],
                     country_list[37:length(country_list)])

for (i in 1:length(country_list)) {

	# By pages
	png(filename=paste0(out_dir, "ad2_by_adm0_", min(year_list), "-", max(year_list), "_pg_", i, ".png"), 
	    type="cairo",
	    units="in", 
	    width=12, 
	    height=7, 
	    pointsize=12, 
	    res=300)

	ad2_page_df <- subset(ad2_df, ADM0_NAME %in% country_list[[i]])

	gg <- ggplot() +
					geom_point(data = ad2_page_df,
					           aes(x = mean_0, xmin = lower_0, xmax = upper_0,
					               y = mean_1, ymin = lower_1, ymax = upper_1,
					               size = pop, color = p_above, shape = sig),
					           alpha = 0.6) +
					geom_abline(alpha = 0.6) +
					geom_abline(slope = 0, intercept = 0.8, 
					            color = "black", alpha = 0.6,
					            linetype = "dotted") +
				#	geom_errorbar(alpha = 0.01) +
				#	geom_errorbarh(alpha = 0.0) +
					scale_size_area() +
				#	scale_alpha_discrete(range = c(0.2, 0.8)) +
				#	scale_color_manual(values = carto_colors) +
				#	scale_colour_distiller(palette = "Reds") +
					scale_color_gradientn(colors = rev(carto_diverging), limits = c(0,1)) +
					facet_wrap(~ADM0_NAME, ncol = 6) +
					theme_minimal() +
					theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
					scale_shape_manual(values = c(1,19)) +
					labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
					     y = paste0("Mean DPT3 Coverage: ", max(year_list)),
					     size = "Number of children \n< 5 years of age",
					     color = paste0("Probability of \nDPT3 coverage > 80% \n(", max(year_list), ")"),
					     shape = "High (>95%) probability \nof true change 2000-2016") +
					scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
					scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
					theme(panel.spacing.x=unit(1, "lines")) +
					coord_equal()

	print(gg)

	dev.off()

}

# Loop over countries

for (adname in unique(ad2_df$ADM0_NAME)) {
	country_df <- subset(ad2_df, ADM0_NAME == adname)

	png(filename=paste0(out_dir, "ad2_by_adm0_", adname, "_", min(year_list), "-", max(year_list), ".png"), 
	    type="cairo",
	    units="in", 
	    width=8, 
	    height=6, 
	    pointsize=12, 
	    res=300)

	gg <- ggplot() +
					geom_point(data = country_df,
					           aes(x = mean_0, xmin = lower_0, xmax = upper_0,
					               y = mean_1, ymin = lower_1, ymax = upper_1,
					               size = pop, color = p_above, shape = sig),
					           alpha = 0.6) +
					geom_abline(alpha = 0.6) +
					geom_abline(slope = 0, intercept = 0.8, 
					            color = "black", alpha = 0.6,
					            linetype = "dotted") +
				#	geom_errorbar(alpha = 0.01) +
				#	geom_errorbarh(alpha = 0.0) +
					scale_size_area() +
				#	scale_alpha_discrete(range = c(0.2, 0.8)) +
				#	scale_color_manual(values = carto_colors) +
				#	scale_colour_distiller(palette = "Reds") +
					scale_color_gradientn(colors = rev(carto_diverging), limits = c(0,1)) +
				#	facet_wrap(~ADM0_NAME, ncol = 6) +
					theme_minimal() +
					theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
					scale_shape_manual(values = c(1,19)) +
					labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
					     y = paste0("Mean DPT3 Coverage: ", max(year_list)),
					     size = "Number of children \n< 5 years of age",
					     color = paste0("Probability of \nDPT3 coverage > 80% \n(", max(year_list), ")"),
					     shape = "High (>95%) probability \nof true change 2000-2015") +
					scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
					scale_y_continuous(expand = c(0, 0), limits = c(0,1))

	print(gg)

	dev.off()

}

#################

png(filename=paste0(out_dir, "ad2_sdi_change", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=300)

gg <- ggplot() +
				geom_point(data = ad2_df,
				           aes(x = mean_0, xmin = lower_0, xmax = upper_0,
				               y = mean_1, ymin = lower_1, ymax = upper_1,
				               size = pop, color = sdi_change, shape = sig),
				           alpha = 0.6) +
			#	geom_errorbar(alpha = 0.01) +
			#	geom_errorbarh(alpha = 0.0) +
				scale_size_area() +
			#	scale_alpha_discrete(range = c(0.2, 0.8)) +
			#	scale_color_manual(values = carto_colors) +
				scale_colour_distiller(palette = "Blues") +
			#	facet_wrap(~ADM0_NAME) +
				theme_classic() +
				scale_shape_manual(values = c(1,19)) +
				labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
				     y = paste0("Mean DPT3 Coverage: ", max(year_list)),
				     size = "Number of children \n< 5 years of age",
				     color = paste0("SDI (", max(year_list), ")"),
				     shape = "High (>95%) probability \nof change 2000-2015") +
				geom_abline() +
				scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
				scale_y_continuous(expand = c(0, 0), limits = c(0,1))

print(gg)

dev.off()

# Try faceted
ad2_df[ADM0_NAME == "Democratic Republic of the Congo", ADM0_NAME := "DRC"]
ad2_df[ADM0_NAME == "Sao Tome and Principe", ADM0_NAME := "STP"]
ad2_df[ADM0_NAME == "United Republic of Tanzania", ADM0_NAME := "Tanzania"]
ad2_df[ADM0_NAME == "Central African Republic", ADM0_NAME := "CAR"]

png(filename=paste0(out_dir, "ad2_sdi_by_adm0_", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=11, 
    height=14, 
    pointsize=12, 
    res=300)

# gg <- ggplot(data = ad2_df, 
#              aes(x = mean_0, xmin = lower_0, xmax = upper_0,
#                  y = mean_1, ymin = lower_1, ymax = upper_1,
#                  size = pop, color = sdi_cut, alpha = sig)) +
# 				geom_point() +
# 			#	geom_errorbar(alpha = 0.01) +
# 			#	geom_errorbarh(alpha = 0.0) +
# 				scale_size_area() +
# 				scale_alpha_discrete(range = c(0.2, 0.8)) +
# 			#	scale_color_manual(values = carto_colors) +
# 				scale_colour_brewer(palette = "Set1") +
# 				facet_wrap(~ADM0_NAME, ncol = 6) +
# 				theme_minimal() +
# 				theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 				# theme(panel.spacing.x=unit(1.3, "lines"),panel.spacing.y=unit(1, "lines")) +
# 				labs(x = paste0("Year: ", min(year_list)),
# 				     y = paste0("Year: ", max(year_list)),
# 				     size = "Number of children \n< 5 years of age",
# 				     color = paste0("SDI (", max(year_list), ")"),
# 				     alpha = "High (>95%) probability \nof change 2000-2015") +
# 				geom_abline() +
# 				coord_equal()

gg <- ggplot() +
		geom_point(data = ad2_df,
		           aes(x = mean_0, xmin = lower_0, xmax = upper_0,
		               y = mean_1, ymin = lower_1, ymax = upper_1,
		               size = pop, color = sdi_cut, shape = sig),
				           alpha = 0.6) +
	#	geom_errorbar(alpha = 0.01) +
	#	geom_errorbarh(alpha = 0.0) +
		scale_size_area() +
	#	scale_alpha_discrete(range = c(0.2, 0.8)) +
	#	scale_color_manual(values = carto_colors) +
		scale_colour_brewer(palette = "Set1") +
	#	facet_wrap(~ADM0_NAME) +
		facet_wrap(~ADM0_NAME, ncol = 6) +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		scale_shape_manual(values = c(1,19)) +
		labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
		     y = paste0("Mean DPT3 Coverage: ", max(year_list)),
		     size = "Number of children \n< 5 years of age",
		     color = paste0("SDI (", max(year_list), ")"),
		     shape = "High (>95%) probability \nof change 2000-2015") +
		geom_abline() +
		scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
		scale_y_continuous(expand = c(0, 0), limits = c(0,1))

print(gg)

dev.off()

# Try by region

for (reg in unique(ad2_df$region)) {

	png(filename=paste0(out_dir, "ad2_by_adm0_", reg, "_", min(year_list), "-", max(year_list), ".png"), 
	    type="cairo",
	    units="in", 
	    width=10, 
	    height=8, 
	    pointsize=12, 
	    res=300)

	gg <- ggplot(data = ad2_df[region == reg], 
	             aes(x = mean_0, xmin = lower_0, xmax = upper_0,
	                 y = mean_1, ymin = lower_1, ymax = upper_1,
	                 size = pop, color = ADM0_NAME, alpha = sig)) +
					geom_point() +
				#	geom_errorbar(alpha = 0.01) +
				#	geom_errorbarh(alpha = 0.0) +
					scale_size_area() +
				#	scale_alpha_discrete(range = c(0.2, 0.8)) +
					scale_color_manual(values = big_colors) +
				#	scale_colour_brewer(palette = "Set1") +
				#	facet_wrap(~ADM0_NAME) +
					theme_bw() +
					labs(x = paste0("Year: ", min(year_list)),
					     y = paste0("Year: ", max(year_list)),
					     size = "Number of children \n< 5 years of age",
					     color = paste0("SDI (", max(year_list), ")"),
					     alpha = "High (>95%) probability \nof change 2000-2015") +
					geom_abline() +
					coord_equal()
}


# Now create plot of mean change vs. mean 2000
ad2_differences <- fread(paste0(in_dir, "pred_derivatives/admin_summaries/",
                   					    indicator, "_admin_2_raked_diff_", 
                   					    min(year_list), "-", max(year_list),".csv"))

ad2_differences <- subset(ad2_differences, select = c("ADM2_CODE", "mean", "upper", "lower"))
setnames(ad2_differences, c("mean", "upper", "lower"), c("mean_diff", "upper_diff", "lower_diff"))

ad2_df <- merge(ad2_df, ad2_differences, by = "ADM2_CODE", all.x = T, all.y = F)

# Plot
png(filename=paste0(out_dir, "ad2_diff_", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=300)

# gg <- ggplot(data = ad2_df, 
#              aes(x = mean_0, xmin = lower_0, xmax = upper_0,
#                  y = mean_diff, ymin = lower_diff, ymax = upper_diff,
#                  size = pop, color = sdi_cut, alpha = sig)) +
# 				geom_point() +
# 			#	geom_errorbar(alpha = 0.01) +
# 			#	geom_errorbarh(alpha = 0.0) +
# 				scale_size_area() +
# 				scale_alpha_discrete(range = c(0.2, 0.8)) +
# 			#	scale_color_manual(values = carto_colors) +
# 				scale_colour_brewer(palette = "Set1") +
# 			#	facet_wrap(~ADM0_NAME) +
# 				theme_classic() +
# 				labs(x = paste0("Year: ", min(year_list)),
# 				     y = paste0("Year: ", max(year_list)),
# 				     size = "Number of children \n< 5 years of age",
# 				     color = paste0("SDI (", max(year_list), ")"),
# 				     alpha = "High (>95%) probability \nof change 2000-2015") +
# 				geom_abline(intercept = 0, slope = 0) +
# 				geom_abline(intercept = 1, slope = -1) +
# 				ylim(NA,1) + xlim(0,1)

gg <- ggplot() +
				geom_point(data = ad2_df,
				           aes(x = mean_0, xmin = lower_0, xmax = upper_0,
				               y = mean_diff, ymin = lower_diff, ymax = upper_diff,
				               size = pop, color = sdi_cut, shape = sig),
				           alpha = 0.6) +
			#	geom_errorbar(alpha = 0.01) +
			#	geom_errorbarh(alpha = 0.0) +
				scale_size_area() +
			#	scale_alpha_discrete(range = c(0.2, 0.8)) +
			#	scale_color_manual(values = carto_colors) +
				scale_colour_brewer(palette = "Set1") +
			#	facet_wrap(~ADM0_NAME) +
				theme_classic() +
				scale_shape_manual(values = c(1,19)) +
				labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
				     y = paste0("Change in DPT3 Coverage: ", min(year_list), " - ", max(year_list)),
				     size = "Number of children \n< 5 years of age",
				     color = paste0("SDI (", max(year_list), ")"),
				     shape = "High (>95%) probability \nof change 2000-2015") +
				geom_abline(intercept = 0, slope = 0) +
				geom_abline(intercept = 1, slope = -1, color = "darkred") +
				scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
				scale_y_continuous(expand = c(0, 0), limits = c(NA,1))

print(gg)

dev.off()

# Plot actual / max possible change
ad2_df[mean_diff > 0, max_change := (1-mean_0)]
ad2_df[mean_diff < 0, max_change := (mean_0)]
ad2_df[,diff_rel_max := mean_diff / max_change]

png(filename=paste0(out_dir, "ad2_potential_diff_", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=300)

# gg <- ggplot(data = ad2_df, 
#              aes(x = mean_0, xmin = lower_0, xmax = upper_0,
#                  y = mean_diff, ymin = lower_diff, ymax = upper_diff,
#                  size = pop, color = sdi_cut, alpha = sig)) +
# 				geom_point() +
# 			#	geom_errorbar(alpha = 0.01) +
# 			#	geom_errorbarh(alpha = 0.0) +
# 				scale_size_area() +
# 				scale_alpha_discrete(range = c(0.2, 0.8)) +
# 			#	scale_color_manual(values = carto_colors) +
# 				scale_colour_brewer(palette = "Set1") +
# 			#	facet_wrap(~ADM0_NAME) +
# 				theme_classic() +
# 				labs(x = paste0("Year: ", min(year_list)),
# 				     y = paste0("Year: ", max(year_list)),
# 				     size = "Number of children \n< 5 years of age",
# 				     color = paste0("SDI (", max(year_list), ")"),
# 				     alpha = "High (>95%) probability \nof change 2000-2015") +
# 				geom_abline(intercept = 0, slope = 0) +
# 				geom_abline(intercept = 1, slope = -1) +
# 				ylim(NA,1) + xlim(0,1)

gg <- ggplot(data = ad2_df,
				           aes(x = mean_0, y = diff_rel_max),
				           alpha = 0.6) +
				geom_point(aes(size = pop, color = sdi_cut, 
				               shape = sig)) +
				geom_abline(intercept = 0, slope = 0) +
				geom_smooth(color = "black", aes(weight = pop)) +
			#	geom_errorbar(alpha = 0.01) +
			#	geom_errorbarh(alpha = 0.0) +
				scale_size_area() +
			#	scale_alpha_discrete(range = c(0.2, 0.8)) +
			#	scale_color_manual(values = carto_colors) +
				scale_colour_brewer(palette = "Set1") +
			#	facet_wrap(~ADM0_NAME) +
				theme_classic() +
				scale_shape_manual(values = c(1,19)) +
				labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
				     y = paste0("Change in DPT3 Coverage \n(relative to maximum possible change)\n", min(year_list), " - ", max(year_list)),
				     size = "Number of children \n< 5 years of age",
				     color = paste0("SDI (", max(year_list), ")"),
				     shape = "High (>95%) probability \nof change 2000-2015") +
				scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
				scale_y_continuous(expand = c(0, 0), limits = c(-1,1)) 
print(gg)

dev.off()

# Another relative difference plot

ad2_df[, rel_diff := mean_diff/(1-mean_0)]

# Plot
png(filename=paste0(out_dir, "ad2_diff_rel_", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=300)

gg <- ggplot(data = ad2_df, 
             aes(x = mean_0, y = rel_diff,
                 size = pop, color = sdi_cut, alpha = sig)) +
				geom_point() +
			#	geom_errorbar(alpha = 0.01) +
			#	geom_errorbarh(alpha = 0.0) +
				scale_size_area() +
				scale_alpha_discrete(range = c(0.2, 0.8)) +
			#	scale_color_manual(values = carto_colors) +
				scale_colour_brewer(palette = "Set1") +
			#	facet_wrap(~ADM0_NAME) +
				theme_classic() +
				labs(x = paste0("Year: ", min(year_list)),
				     y = paste0("Relative improvement"),
				     size = "Number of children \n< 5 years of age",
				     color = paste0("SDI (", max(year_list), ")"),
				     alpha = "High (>95%) probability \nof change 2000-2015") +
				geom_abline(intercept = 0, slope = 0) 

print(gg)

dev.off()