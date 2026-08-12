library(dplyr)
library(tidyr)

tar_load(regen_metrics)
tar_load(metric3d_all)

# metric3d_all_plot=metric3d_all %>% 
#   filter(subplot=="plot") %>% 
#   pivot_longer(cols =  c(matches("h_"),"coef_lin")) 
# 
# metric3d_all %>% 
#   filter(subplot!="plot") %>% 
#   pivot_longer(cols = c(matches("h_"),"coef_lin")) %>% 
#   ggplot(aes(slice_num,value,color=plot))+
#   geom_line(aes(group=interaction(plot,subplot)),alpha = 0.4)+
#   geom_point(data=metric3d_all_plot,aes(slice_num,value,color=plot))+
#   geom_line(data=metric3d_all_plot,aes(slice_num,value,color=plot,group=interaction(plot,subplot)))+
#   facet_grid(name~plot,scale="free_y")+
#   theme_bw()

# 
# reg_ext=regen_metrics_subplot %>% 
#   select(plot,subplot,richness,abundance,H,browsing) %>% 
#   pivot_longer(cols=-c(plot,subplot),names_to = "y",values_to = "y_val")
# 
# 
# metric3d_all %>% 
#   select(plot,subplot,slice_num,rb_sd) %>% filter(slice_num==1) %>% 
#   filter(subplot!="plot") %>%
#   left_join(reg_ext,by=c("plot","subplot")) %>% 
#   distinct() %>%  
#   ggplot(aes(rb_sd,y_val,col=plot))+
#   geom_point()+geom_smooth(aes(group=plot))+
#   facet_wrap(~y,scale="free_y",ncol=4)

# 
# regen_metrics %>% 
#   rename(plot=Plot,
#          subplot=Subplot) %>% 
#   ggplot(aes(browsing,abundance,color=plot))+
#   geom_point()
# 
# 
# # variable n_seedling is useless, because seedling are the most abundant class
# regen_metrics %>% 
#   ggplot(aes(abundance,n_seedling))+
#   geom_point()
# 
# regen_metrics_plot=regen_metrics %>% 
#   select(-Species) %>% 
#   group_by(Plot,Subplot) %>% 
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)),
#             .groups = "drop")  


# effect of browsing on structural variables
# regen_metrics_plot %>% 
#   pivot_longer(cols=-c("Plot","Subplot","browsing")) %>% 
#   mutate(name=factor(name,levels = c("abundance","abundance_cv","n_seedling","richness","richness_sd",
#                                      "Height_increment_mn_hclass2","Height_increment_mn_hclass3",
#                                      "Height_increment_cv_hclass2","Height_increment_cv_hclass3",
#                                      "H",
#                                      "H_hclass1","H_hclass2","H_hclass3","H_hclass4","H_hclass5"))) %>% 
#   ggplot(aes(browsing,value))+
#   geom_point()+
#   geom_smooth(method="gam")+
#   facet_wrap(~name,scale="free_y",ncol=5)+
#   theme_bw()

# analysis per species
browsing_sp=tar_read(regen_df) %>%
  filter(Hclass > 1) %>%
  group_by(Plot, Subplot,Species) %>%
  summarise(browsing_sp = sum(Browsing == "Y") / n(), .groups = "drop")
browsing_plot=tar_read(regen_df) %>%
  filter(Hclass > 1) %>%
  group_by(Plot, Subplot) %>%
  summarise(browsing = sum(Browsing == "Y") / n(), .groups = "drop")
browsing_plot %>% 
  left_join(browsing_sp) %>% 
  mutate(browsing_dif=browsing-browsing_sp) %>% 
  left_join(regen_metrics %>% 
              select(Plot,Subplot,Species,matches("Height_"))) %>% 
  pivot_longer(cols=-c("Plot","Subplot","Species",matches("browsing"))) %>% 
  ggplot(aes(browsing_dif,value,color=Species))+
  geom_point()+
  geom_smooth(aes(group=NA),method="gam")+
  facet_wrap(~name,scale="free_y")+
  theme_bw()

# comparing 3d struct and light
metric3d_all_plot %>% 
  ggplot(aes(h_sd,value))+
  geom_point()+geom_line()+
  facet_grid(name~slice_num,scale="free")



## explore range of variables
regen_metrics %>% 
  select(Plot,Subplot,richness,abundance,H,Height_increment_mn_hclass2) %>% 
  rename(plot=Plot,
         subplot=Subplot) %>% 
  group_by(plot,subplot) %>% 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>% 
  mutate(subplot=paste0("subplot_",subplot)) %>% 
  pivot_longer(cols=-c("plot","subplot")) %>% 
  mutate(slice_num=NA) %>% 
  filter(plot%in%metric3d_all$plot)->regen_ext
metric3d_all %>% 
  select(plot,subplot,slice_num,rb_mean,rb_sd,rb_cv) %>% 
  filter(slice_num%in%c(1,4)) %>% 
  pivot_longer(cols=c("rb_mean","rb_sd","rb_cv")) %>% 
  bind_rows(regen_ext) %>% 
  mutate(name=factor(name,level=c("rb_mean","rb_sd","rb_cv",
                                  "abundance","richness","H","Height_increment_mn_hclass2"))) %>% 
  ggplot(aes(value,plot))+
  geom_boxplot(aes(color=as.factor(slice_num)))+
  facet_wrap(~name,scale="free_x",ncol=3)+
  theme_bw()+
  labs(color="Height slice",x="Radiative budget metric or Structural metric",y="")



metric3d_all %>% 
  select(plot,subplot,slice_num,rb_mean,rb_sd,rb_cv) %>% 
  filter(slice_num%in%c(1,4)) %>% 
  left_join(regen_ext %>% 
              pivot_wider(names_from = name,values_from = value) %>% 
              select(-slice_num),by=c("plot","subplot")) %>% 
  ggplot(aes(rb_mean,richness))+
  geom_point(aes(color=plot))+
  geom_smooth(method="lm")+
  facet_wrap(~slice_num,scale="free")


# evolution of radiative budget throughout the day and the plot/subplots
summary_plot=tar_read(metric3d_df_GER11)

summary_plot %>% 
  ggplot(aes(time,rb_mean,color=as.factor(slice_num)))+
  geom_point()+
  geom_line()+
  facet_wrap(~subplot)+
  scale_x_continuous(breaks=seq(8, 18, 1))




# check the light distribution for subplot 1
pc_norm=loadRData(tar_read(pc_norm_GER11))
list_rb=loadRData(tar_read(rb_norm_GER11))
bbox=tar_read(bbox_GER11)
subplot_extent=tar_read(subplot_extent_GER11)
plot_name="GER11"


rb_norm=list_rb[[3]]
pc_slice_chm <- rasterize_canopy(pc_norm, res = 0.1, algorithm = p2r())
slice_list = seq(from = 0.5, by = 1, length.out = 5)
list_g=vector("list",length=4)
list_s=vector("list",length=4)
for(i in 1:4){
  slice_low=slice_list[i]
  slice_up=slice_list[i+1]
  rb_slice=rast(apply(rb_norm[,,(slice_low/0.25):(slice_up/0.25)], c(1, 2), sum))
  ext(rb_slice)  <- ext(pc_slice_chm)
  crs(rb_slice)  <- crs(pc_slice_chm)
  rb_slice_matched <- resample(rb_slice, pc_slice_chm, method = "bilinear")
  names(rb_slice_matched)="rb"
  rb_slice_matched[rb_slice_matched<quantile(values(rb_slice_matched),
                                             na.rm=TRUE,probs=0.005)] <-NA
  plot_ext=subplot_extent$subplot_1
  plot_ext_trans <- shift(vect(plot_ext), 
                          dx = -bbox$xmin, 
                          dy = -bbox$ymin)
  pc_sub_slice_chm=crop(pc_slice_chm,plot_ext_trans)
  rb_sub_slice_matched=crop(rb_slice_matched,plot_ext_trans)
  
  list_g[[i]]=ggplot() +
    geom_spatraster(data = rb_slice_matched, aes(fill = rb)) 
  list_s[[i]]=ggplot() +
    geom_spatraster(data = rb_sub_slice_matched, aes(fill = rb))
}
