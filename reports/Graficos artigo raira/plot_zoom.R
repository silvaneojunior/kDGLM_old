plot_zoom=function(p1,mini_range,mini_place){
  p2=p1+
    coord_cartesian(ylim=c(mini_range$y[1],mini_range$y[2]),xlim=c(mini_range$x[1],mini_range$x[2]))+
    guides(color=FALSE,shape=FALSE,linetype=FALSE,fill=FALSE)+
    scale_x_continuous('',expand=c(0,0))+
    scale_y_continuous('',expand=c(0,0))+
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank(),
          plot.margin=unit(rep(0,4), "cm"))

  data_mini_place=data.frame('xmin'=mini_place$x[1],'xmax'=mini_place$x[2],
                             'ymin'=mini_place$y[1],'ymax'=mini_place$y[2])
  data_mini_range=data.frame('xmin'=mini_range$x[1],'xmax'=mini_range$x[2],
                             'ymin'=mini_range$y[1],'ymax'=mini_range$y[2])

  p3=p1+
    geom_path(data=data.frame(x = c(mini_range$x[1],
                                    mini_place$x[1],
                                    mini_range$x[2],
                                    mini_place$x[2]),
                              y=c(mini_range$y[1],
                                  mini_place$y[1],
                                  mini_range$y[1],
                                  mini_place$y[1])
                              ,grp=c(1,1,2,2)),
              aes(x,y,group=grp),
              linetype='dashed') +
    geom_path(data=data.frame(x = c(mini_range$x[1],
                                    mini_place$x[1],
                                    mini_range$x[2],
                                    mini_place$x[2]),
                              y=c(mini_range$y[2],
                                  mini_place$y[2],
                                  mini_range$y[2],
                                  mini_place$y[2])
                              ,grp=c(1,1,2,2)),
              aes(x,y,group=grp),
              linetype='dashed') +
    geom_rect(data=data_mini_place,
              aes(xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax), color='black', fill='white', linetype='solid', alpha=1,size=1.5) +
    annotation_custom(ggplotGrob(p2),
                      xmin = mini_place$x[1],
                      xmax = mini_place$x[2],
                      ymin = mini_place$y[1],
                      ymax = mini_place$y[2]) +
    geom_rect(data=data_mini_range,
              aes(xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax), color='black', linetype='solid', alpha=0)
  return(p3)
}

