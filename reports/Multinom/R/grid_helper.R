# library(tidyr)
# library(hms)
# library(beepr)

guys=c(1,2,3,4,5)
ranges=list(c(1:10),c(1:10),c(1:10),c(1:10),c(1:10))

run_grid=function(ranges,func_eval){
  n_var=length(ranges)
  sizes=c()
  grid=list()
  for(i in c(1:n_var)){
    sizes=c(sizes,length(ranges[[i]]))
    grid[[i]]=ranges[[i]]
  }
  total=prod(sizes)
  ref_data=as.data.frame(matrix(0,prod(sizes),n_var+1))

  set_values=function(index,value){
    ref_data[index,]<<-c(value,func_eval(value))
  }
  count=0
  init=Sys.time()
  while(count<=total){
    perc=count/total
    cur_time=Sys.time()

    qtd1=min(49,floor(49*perc))
    qtd2=49-qtd1

    cat(paste0('[',paste0(rep('=',qtd1),collapse = ''),'>',paste0(rep(' ',qtd2),collapse = ''),']  ',
               paste(formatC(count,format='d',big.mark=','),'/',formatC(total,format='d',big.mark=',')) %>% paste0(' (',(100*perc) %>% round(2) %>% format(nsmall=2) %>% paste0('%)')),
               ' - ETA: ',((1-perc)*difftime(cur_time, init, unit="secs")/perc) %>% as.numeric %>% round %>% hms,
               '\r'))

    marker=count
    value=c(1:n_var)
    for(i in c(1:n_var)){
      value[i]=grid[[i]][marker%%sizes[i]+1]
      marker=marker%/%sizes[i]
    }
    set_values(count,value)
    count=count+1
  }

  beep()
  return(ref_data)
}

#run_grid(ranges,function(x){rnorm(1)})
