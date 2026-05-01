library(dplyr)
library(tidyr)

tar_load(regen_metrics)
tar_load(metric3d_all)

metric3d_all_plot=metric3d_all %>% 
  filter(subplot=="plot") %>% 
  pivot_longer(cols =  c(matches("h_"),"coef_lin")) 

metric3d_all %>% 
  filter(subplot!="plot") %>% 
  pivot_longer(cols = c(matches("h_"),"coef_lin")) %>% 
  ggplot(aes(slice_num,value,color=plot))+
  geom_line(aes(group=interaction(plot,subplot)),alpha = 0.4)+
  geom_point(data=metric3d_all_plot,aes(slice_num,value,color=plot))+
  geom_line(data=metric3d_all_plot,aes(slice_num,value,color=plot,group=interaction(plot,subplot)))+
  facet_grid(name~plot,scale="free_y")+
  theme_bw()


reg_ext=regen_metrics %>% 
  select(Plot,Subplot,richness,abundance,H,browsing) %>% 
  rename(plot=Plot,
         subplot=Subplot) %>% 
  mutate(subplot=paste0("subplot_",subplot)) %>% 
  pivot_longer(cols=-c(plot,subplot),names_to = "y",values_to = "y_val")


metric3d_all %>% 
  select(plot,subplot,slice_num,rb_sd) %>% filter(slice_num==1) %>% 
  filter(subplot!="plot") %>%
  left_join(reg_ext,by=c("plot","subplot")) %>% 
  distinct() %>%  
  ggplot(aes(rb_sd,y_val,col=plot))+
  geom_point()+geom_smooth(aes(group=plot))+
  facet_wrap(~y,scale="free_y",ncol=4)


regen_metrics %>% 
  rename(plot=Plot,
         subplot=Subplot) %>% 
  ggplot(aes(browsing,abundance,color=plot))+
  geom_point()


# variable n_seedling is useless, because seedling are the most abundant class
regen_metrics %>% 
  ggplot(aes(abundance,n_seedling))+
  geom_point()

regen_metrics_plot=regen_metrics %>% 
  select(-Species) %>% 
  group_by(Plot,Subplot) %>% 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")  


# effect of browsing on structural variables
regen_metrics_plot %>% 
  pivot_longer(cols=-c("Plot","Subplot","browsing")) %>% 
  mutate(name=factor(name,levels = c("abundance","abundance_cv","n_seedling","richness","richness_sd",
                                     "Height_increment_mn_hclass2","Height_increment_mn_hclass3",
                                     "Height_increment_cv_hclass2","Height_increment_cv_hclass3",
                                     "H",
                                     "H_hclass1","H_hclass2","H_hclass3","H_hclass4","H_hclass5"))) %>% 
  ggplot(aes(browsing,value))+
  geom_point()+
  geom_smooth(method="gam")+
  facet_wrap(~name,scale="free_y",ncol=5)+
  theme_bw()

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
