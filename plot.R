library(dplyr)
library(tidyr)

tar_load(regen_metrics)
tar_load(metric3d_all)

metric3d_all %>% 
  pivot_longer(cols = matches("rb_")) %>% 
  ggplot(aes(slice_num,value,color=plot))+
  geom_point()+geom_line()+
  facet_grid(name~subplot,scale="free_y")


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
