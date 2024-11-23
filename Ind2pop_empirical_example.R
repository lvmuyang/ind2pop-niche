# Working example for individual to population niche
# Libraries
library(tidyverse)
library(lubridate)
library(sf)
library(glue)
library(lme4)
library(ggplot2)
library(viridisLite)
library(moments)
library(ctmm)
library(raster)
library(assertthat)
library(RColorBrewer)
library(rnaturalearth)
library(ggmap)
library(patchwork)
library(ggpubr)
library(acid)

# Custom Functions
source("./Ind2pop_function.R")

#----   Load Data   ----#
message("Loading data...")
sp_anno <- read.csv(file = "gadwall_annotated.csv")
sp_name <- "Gadwall"

sp_anno <- read.csv(file = "elephants_annotated.csv")
sp_name <- "Elephant"
# 2008.10.30 - 2014.03.27

# Environmental Raster
if(sp_name == "Gadwall") lst <- raster("./modis_lst_mean_2009-2010.tif") 
if(sp_name == "Elephant") lst <- raster("./modis_lst_mean_2008-2014.tif") 
values(lst) <- values(lst)*.02-273.15 # convert to C
# convert K to C
# sp_anno <- sp_anno %>% 
#   mutate(value = value*.02-273.15)
warm = 2.4

#----   Perform Analyses - estimate home range with CTMM package   ----#
# This script estimates home ranges for each individual which are later used to 
# set the spatial extend for climate vuinerability analyses.
inds <- unique(sp_anno$individual.local.identifier)

# initialize list to collect results
out <- list()
n_max <- 500
for(i in 1:length(inds)){
  set.seed(123)
  message(glue("Starting ctmm process for Individual {inds[i]}..."))
  #format data as `telemetry` object for ctmm
  evt_tmp <- sp_anno %>% 
    filter(individual.local.identifier == inds[i])
  if(nrow(evt_tmp)>=n_max) evt_tmp <- evt_tmp[sort(sample(1:nrow(evt_tmp),n_max)),]
  
  if(nrow(evt_tmp) == 0){
    message("Insufficient records during season, moving to next...")
  } else{
    message(glue("{nrow(evt_tmp)} records found..."))
    evt_telem <- as.telemetry(evt_tmp)
  } 
  
  # get initial acf guess
  guess <- ctmm.guess(evt_telem, interactive = F)
  
  # acf model selection using guess as init
  fit <- ctmm.select(evt_telem, guess)
  
  # fit 
  akde <- akde(evt_telem, fit)
  
  message(glue("Gathering output for Individual {inds[i]}..."))
  # write the akde object to list
  tmp_out <- list()
  tmp_out[["Individual"]] <- as.character(inds[i])
  tmp_out[["aKDE Object"]] <- akde
  tmp_out[["ctmm Fit"]] <- fit
  tmp_out[["PDF Raster"]] <- raster(akde, DF = "PDF")
  tmp_out[["SF Polygon"]] <- as.sf(akde, level.UD=.95)
  tmp_out[["Clipped Raster"]] <- terra::mask(tmp_out[["PDF Raster"]], tmp_out[["SF Polygon"]])
  
  out[[i]] <- tmp_out
  
  message(glue("Individual {inds[i]} complete..."))
}

message("Saving aKDEs to out/akde_out.rdata...")
save(out, file = paste0("out/",sp_name,"_akde_out.rdata"))

#----   Perform Analyses - climate change vulnerability   ----#
# This script estimates individual climate change vulnerabilities
# aKDEs
#load("out/akde_out.rdata")
akdes <- out

#----   Perform Analyses   ----#
inds <- unique(sp_anno$individual.local.identifier)

# initialize list to collect results
out <- list()
#i <- 1
pdf(paste0(sp_name,"_individual_niches.pdf"),height = 4,width = 12)
for(i in 1:length(inds)){
  
  message(glue("Calculating vulnerability for Individual {inds[i]}..."))
  
  dat_ind <- sp_anno %>% 
    filter(individual.local.identifier == inds[i])
  #if(i==4) dat_ind <- dat_ind[-c(1:91),]
  begin_date <- dat_ind$date[1]
  end_date <- dat_ind$date[length(dat_ind$date)]
  
  
  akde <- akdes[[i]]
  hr_poly <- akde$`SF Polygon`[glue("{inds[i]} 95% est"),] %>% 
    st_transform(crs = projection(lst))
  
  # Check that we got the right ind...
  assert_that(akde[["Individual"]] == inds[i])
  
  # Make LST raster for individual home range
  ind_lst <- lst %>% terra::crop(hr_poly)
  #ind_lst <- lst %>% terra::mask(hr_poly)
  vals <- values(ind_lst)
  #vals <- values(lst)
  
  coordinates(dat_ind) <- ~lng + lat
  lst_pj <- projectRaster(lst,akde$`PDF Raster`)
  vals_ind_lst <- values(lst_pj)
  weights_ind_lst <- values(akde$`PDF Raster`)/sum(values(akde$`PDF Raster`))
  dat_mods_ind <- cbind.data.frame(vals_ind_lst,weights_ind_lst)
  dat_mods_ind <- dat_mods_ind[!is.na(dat_mods_ind$vals_ind_lst),]
  dat_mods_ind$weights_ind_lst <- dat_mods_ind$weights_ind_lst/sum(dat_mods_ind$weights_ind_lst)
  
  d_ind_lst <- density(x = dat_mods_ind$vals_ind_lst, weights = dat_mods_ind$weights_ind_lst,
                       from = 15,
                       to = 45,
                       bw = 0.5)
  d_ind_lst_track <- density(x = na.omit(dat_ind$value), weights = NULL,
                             from = 15,
                             to = 45,
                             bw = 0.5)
  d_lst <- density(x = na.omit(vals_ind_lst), weights = NULL,
                   from = 15,
                   to = 45,
                   bw = 0.5)
  par(mfrow=c(1,3))
  plot(lst,main = inds[i])
  plot(hr_poly,add=T)
  plot(ind_lst,main = paste0(begin_date," to ",end_date))
  plot(hr_poly,add=T)
  plot(dat_ind,add=T,col = 2)
  plot(d_ind_lst$y~d_ind_lst$x,type='l',ylab = "Density", xlab = "LST",main = inds[i])
  lines(d_lst$y~d_lst$x,col=2)
  lines(d_ind_lst_track$y~d_ind_lst_track$x,col=4)
  
  f_ind_lst <- d_ind_lst$y/d_lst$y
  f_ind_lst[is.infinite(f_ind_lst)] <- 0
  #f_ind_lst <- f_ind_lst/sum(f_ind_lst,na.rm = T)
  f_ind_lst <- f_ind_lst/max(f_ind_lst,na.rm = T)/3
  lines(f_ind_lst~d_ind_lst_track$x,col=5)
  legend("topright",c("CTMM-derived","LST in home range","LST tracked"),lty=1,col=c(1,2,4),cex=1)
  
  
  # plot(f_ind_lst~d_lst$x,col=4,type ='l')
  # approx <- approxfun(d$x, d$y)
  # Look up the probability associated with some value x.
  # prob <- approx(x)
  
  cur <- sum(sapply(vals, getUDVal_wt, dat = dat_mods_ind$vals_ind_lst, w = dat_mods_ind$weights_ind_lst), na.rm = T)
  fut <- sum(sapply(vals+warm, getUDVal_wt, dat = dat_mods_ind$vals_ind_lst, w = dat_mods_ind$weights_ind_lst), na.rm = T)
  
  # cur <- sum(sapply(vals, getUDVal, dat = dat_ind$value), na.rm = T)
  # fut <- sum(sapply(vals+warm, getUDVal, dat = dat_ind$value), na.rm = T)
  
  w <- fut/cur
  
  out[[i]] <- data.frame(individual = inds[i],
                         fut_weight = w,
                         temp = d_ind_lst$x,
                         suit = d_ind_lst$y)
  
}  
dev.off()
df_out <- do.call("rbind", out)

saveRDS(df_out, paste0("out/",sp_name,"_fut-weights.rds"))

#----   Perform Analyses   ----#
#-- Init
fig_w <- 3.2
fig_h <- 3.2
if(sp_name == 'Gadwall') plot_temp_range <- c(20,32)
if(sp_name == 'Elephant') plot_temp_range <- c(30,42)


df_out <- filter(df_out,individual!="91891B")
inds_keep <- unique(df_out$individual)

(sp_dens <- ggplot()+
    geom_line(data = df_out, aes(x = temp,
                                 y = suit,
                                 group = individual,
                                 color = fut_weight),
              alpha = 0.7,
    ) +
    # scale_color_gradient2(low = "red", mid = "white", high = "blue")+
    scale_color_viridis_c(direction = 1) +
    xlab("")+
    #xlim(20,32)+
    xlim(plot_temp_range)+
    ylab("")+
    theme_classic() +
    labs(color = "")+
    theme(legend.position = c(.9,.8)))

# Save plot
ggsave(filename = paste0("out/figs/",sp_name,"_dens.png"), sp_dens, width = fig_w, 
       height = fig_h, units = "in")

# future and current niches
current_niche <- 0
future_niche <- 0
for (i in 1:length(inds_keep)) {
  current_niche <- filter(df_out,individual==inds_keep[i])$suit + current_niche
  future_niche <- filter(df_out,individual==inds_keep[i])$suit*
    filter(df_out,individual==inds_keep[i])$fut_weight + future_niche
}
future_niche <- future_niche/sum(unique(df_out$fut_weight))
current_niche <- current_niche/length(unique(df_out$fut_weight))
pop_niche <- cbind.data.frame(suit = c(current_niche,future_niche),
                              temp = filter(df_out,individual %in% inds_keep[1:2])$temp,
                              time = rep(c("current","future"),each = length(current_niche)))

(pop_dens <- ggplot()+
    geom_line(data = pop_niche, aes(x = temp,
                                    y = suit,
                                    #group = time,
                                    color = time)) +
    # scale_color_gradient2(low = "red", mid = "white", high = "blue")+
    #scale_color_viridis_c(direction = 1) +
    xlab("")+
    xlim(plot_temp_range)+
    ylab("")+
    labs(color = "")+
    theme_classic() +
    theme(legend.position = c(0.9,0.8)))

# Save plot
ggsave(filename = paste0("out/figs/",sp_name,"_pop_dens.png"), pop_dens, width = fig_w, 
       height = fig_h, units = "in")


weighted.moments(filter(df_out,individual==inds_keep[i])$temp,current_niche)
weighted.moments(filter(df_out,individual==inds_keep[i])$temp,future_niche)
#---- Bivariate Plot ----#
# for paritioning individual contribution to niche breadth and skewness
message("Creating bivariate plot...")

head(df_out)
ind_para <- matrix(NA , nrow = length(inds_keep), ncol = 5)
colnames(ind_para) <- c("individual","fut_weight","mu","sigma","skew")
ind_para <- as.data.frame(ind_para)
for (i in 1:length(inds_keep)) {
  ind_dat <- filter(df_out, individual == inds_keep[i])
  #skewness(sample(ind_dat$temp,size = 10000000,replace = T,prob = ind_dat$suit))
  #weighted.moments(ind_dat$temp,ind_dat$suit)$w.skew.Stata
  ind_para[i,"individual"] <- inds_keep[i]
  ind_para[i,"fut_weight"] <- unique(ind_dat$fut_weight)
  ind_para[i,"mu"] <- weighted.moments(ind_dat$temp,ind_dat$suit)$weighted.mean
  ind_para[i,"sigma"] <- weighted.moments(ind_dat$temp,ind_dat$suit)$weighted.sd
  ind_para[i,"skew"] <- weighted.moments(ind_dat$temp,ind_dat$suit)$w.skew.Stata
}

ind_conti <- cbind.data.frame(ind_para,individual_contribution(ind_para)) 

(sp_bivar <- ggplot(ind_conti) +
    geom_point(aes(x=marginality_sigma2, y = specialization_sigma2,
                   color = fut_weight), alpha = 1, size = 3)+
    geom_abline(slope = 1, intercept = 0,linetype="dashed",color="grey")+
    scale_size_continuous(name = "Skewness") +
    # scale_color_gradient2(low = "red", mid = "white", high = "blue")+
    scale_color_viridis_c(direction = 1, name = "") +
    # xlim(0,max(ind_conti$marginality_sigma2,ind_conti$specialization_sigma2))+
    # ylim(0,max(ind_conti$marginality_sigma2,ind_conti$specialization_sigma2))+
    xlab("") + 
    ylab("") +
    theme_linedraw()+
    theme_classic() +
    guides(color = guide_colorbar(title.position = "top",
                                  title.hjust = 0.5))+
    theme(legend.position="none", legend.justification = "center",
          legend.text=element_text(size=6),
          legend.title=element_text(size=8))
)

(sp_skew <- ggplot(ind_conti) +
    geom_point(aes(x=marginality_skew, y = specialization_skew,
                   color = fut_weight), alpha = 1, size = 3)+
    geom_vline(xintercept = 0,linetype="dashed",color="grey")+
    geom_hline(yintercept = 0,linetype="dashed",color="grey")+
    scale_size_continuous(name = "Skewness") +
    # scale_color_gradient2(low = "red", mid = "white", high = "blue")+
    scale_color_viridis_c(direction = 1, name = "") +
    xlab("") + 
    ylab("") +
    theme_linedraw()+
    theme_classic() +
    guides(color = guide_colorbar(title.position = "top",
                                  title.hjust = 0.5))+
    theme(legend.position="none", legend.justification = "center",
          legend.text=element_text(size=6),
          legend.title=element_text(size=8))
)
# Save plot
ggsave(filename = paste0("out/figs/",sp_name,"_ind_var.png"), sp_bivar, width = fig_w, 
       height = fig_h, units = "in")
ggsave(filename = paste0("out/figs/",sp_name,"_ind_skew.png"), sp_skew, width = fig_w, 
       height = fig_h, units = "in")

#---- Quick Maps

message("Creating maps...")

world <- ne_countries(returnclass = "sf")

europe <- world[world$continent == "Europe",]
africa <- world[world$continent == "Africa",]
if(sp_name == 'Gadwall'){
  region <- europe
  inset_box <- c(-20,45,30,73)
  inset_range <- c(13.7,14.7,47.5,48.5)
  dens_inset_range <- c(19.5,23,0.5,0.8)
  boxes<-data.frame(maxlat = 50,minlat = 47,maxlong = 15,minlong = 10, id="1")
} 
if(sp_name == 'Elephant'){
  region <- africa
  inset_box <- extent(region)
  inset_range <- c(14,15,-20.5,-19.5)
  dens_inset_range <- c(30,34,0.3,0.5)
  boxes<-data.frame(maxlat = -17,minlat = -21,maxlong = 19,minlong = 12, id="1")
} 
boxes<-transform(boxes, laby=(maxlat +minlat )/2, labx=(maxlong+minlong )/2)

#sp <- read.csv("data/spwall_annotated.csv") %>% 
# sp <- read.csv("spwall_annotated.csv") %>% 
#   left_join(tot2, by = c("individual.local.identifier"="ID")) %>% 
#   st_as_sf(coords = c("lng", "lat"), crs = 4326)
sp <- sp_anno %>%
  #filter(., individual.local.identifier %in% inds_keep) %>%
  left_join(ind_conti, by = c("individual.local.identifier"="individual")) %>% 
  st_as_sf(coords = c("lng", "lat"), crs = 4326)

sp_cent_sf <- st_as_sfc(st_bbox(sp)) %>% 
  st_centroid() 

sp_cent <- sp_cent_sf %>% 
  st_coordinates()

register_google(key = "AIzaSyBhtt2B4pG7bxt69ntHnorlXsBJ0elDz9k")
sp_bg <- get_googlemap(center = sp_cent, zoom = 8, maptype = "terrain", 
                        color = "bw")


sp_inset <- ggplot() +
  geom_sf(data = region)+
  #geom_sf(data = sp_cent_sf, color = "firebrick2", size = 3)+
  geom_rect(data=boxes, aes(xmin=minlong , xmax=maxlong, ymin=minlat, ymax=maxlat ), color="red", fill="transparent") + 
  coord_sf(crs = 4326,
           xlim = inset_box[1:2], 
           ylim = inset_box[3:4])+
  theme_bw(base_size = 12) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.ticks.length = unit(0, "mm"),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(color = "transparent", fill = "transparent"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)
        # panel.spacing = unit(c(0, 0, 0, 0), "null")
  )

(sp_main <- ggmap(sp_bg, darken = c(0.8, "white")) +
    geom_sf(data = sp, aes(color = fut_weight), 
            alpha = 1, size = 1, inherit.aes = F)+
    scale_color_viridis_c(direction = 1) +
    theme_bw(base_size = 12) +
    coord_sf(crs = 4326,
             # xlim = c(10.7, 13),
             # ylim = c(48, 49)
    ) +
    theme(legend.position = "none",
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank(),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))+
    annotation_custom(grob = ggplotGrob(sp_inset),
                      xmin = inset_range[1], xmax = inset_range[2],
                      ymin = inset_range[3], ymax = inset_range[4])
)

# sp_main + inset_element(sp_inset, 0.6, -0.56, 0.99, 1)
ggsave(filename = paste0("out/figs/",sp_name,"_map_main.png"), sp_main, height = fig_h, width = fig_w)
ggsave(filename = paste0("out/figs/",sp_name,"_map_inset.png"), sp_inset, height = 1.25, width = 1.25)



(sp_dens_inset <- sp_dens +
    annotation_custom(grob = ggplotGrob(sp_inset),
                      xmin = dens_inset_range[1], xmax = dens_inset_range[2],
                      ymin = dens_inset_range[3], ymax = dens_inset_range[4]))

ggsave(filename = paste0("out/figs/",sp_name,"_dens_inset.png"), sp_dens_inset, height = fig_h, width = fig_w)

# make raster maps
lst_df <-
  as.data.frame(lst, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit()
colnames(lst_df)[3] <- "mean LST"

# cur <- sum(sapply(vals, getUDVal_wt, dat = dat_mods_ind$vals_ind_lst, w = dat_mods_ind$weights_ind_lst), na.rm = T)
# fut <- sum(sapply(vals+warm, getUDVal_wt, dat = dat_mods_ind$vals_ind_lst, w = dat_mods_ind$weights_ind_lst), na.rm = T)

lst_df$current <- sapply(lst_df$`mean LST`, getUDVal_wt, dat = pop_niche$temp[pop_niche$time == 'current'], w = pop_niche$suit[pop_niche$time == 'current']/sum(pop_niche$suit[pop_niche$time == 'current']))
lst_df$future <- sapply(lst_df$`mean LST` + warm, getUDVal_wt, dat = pop_niche$temp[pop_niche$time == 'future'], w = pop_niche$suit[pop_niche$time == 'future']/sum(pop_niche$suit[pop_niche$time == 'future']))

(
  lst_map <- ggplot(data = lst_df) +
    geom_raster(aes(x = x, y = y, fill = `mean LST`)) +
    scale_fill_viridis_c() +
    theme_void() +
    theme(
      legend.position = "bottom"
    )
)
(
  lst_map_current <- ggplot(data = lst_df) +
    geom_raster(aes(x = x, y = y, fill = `current`)) +
    scale_fill_viridis_c() +
    theme_void() +
    theme(
      legend.position = "bottom"
    )
)
(
  lst_map_future <- ggplot(data = lst_df) +
    geom_raster(aes(x = x, y = y, fill = `future`)) +
    scale_fill_viridis_c() +
    theme_void() +
    theme(
      legend.position = "bottom"
    )
)

ggsave(filename = paste0("out/figs/",sp_name,"_mean_lst.png"), lst_map, height = fig_h, width = fig_w)
ggsave(filename = paste0("out/figs/",sp_name,"_current_pop_suit.png"), lst_map_current, height = fig_h, width = fig_w)
ggsave(filename = paste0("out/figs/",sp_name,"_future_pop_suit.png"), lst_map_future, height = fig_h, width = fig_w)
