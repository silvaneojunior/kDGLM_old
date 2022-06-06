# Shiny stuff
library(shiny)
library(shinyWidgets)
library(shinyjs)
library(shinyjqui)

# Plotting stuff
library(ggplot2)
library(plotly)

# Manipulation  stuff
library(Matrix)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(abind)

# Misc stuff
library(hms)
library(beepr)
library(latex2exp)

# Math stuff
library(rootSolve)
library(MASS)

source('R/main.R')
source('R/ui_helper.R')

final_ui=fluidPage(
  'Generialized Dinamic Linear Models for Exponential family auxiliar aplication',
  useShinyjs(),
  withMathJax(),
  sidebarLayout(
    uiOutput('side_panel'),
    uiOutput('main_panel')
  ))


abind_aux_array=function(...){
  abind(...,along=1)
}

server=function(input, output) {

  modal_prompt=modalDialog(
    #numericInput('column_index', 'Select the index of the column to insert:', value = 1,min=1,max=,step=1),
    fluidRow(
      column(8,
             uiOutput('modal_content')),
      column(4,
             fluidRow(selectInput('kernel','Target variable distribuition',choices=list('Poisson (univariada)','Multinomial','Normal'))),
             fluidRow(materialSwitch('by_column','Variables on rows:')),
             fluidRow(materialSwitch('column_name','Column names:',value=TRUE)),
             fluidRow(materialSwitch('row_name','Row names:',value=TRUE)))
      ),
    footer = actionButton('choosed_column','Select')
  )

  output$modal_content=renderUI({
    selectInput('column_index','Select the column to insert:',
                 choices=names(raw_data()),
                 multiple=kernel_list[[input$kernel]][['multi_var']])
  })

  observeEvent(c(input$y_resp),{
    file=raw_data()
    if(length(dim(file))>1){
      if(dim(file)[2]>1){
        showModal(modal_prompt)
      }
    }
    })

  raw_data=eventReactive(c(input$y_resp,input$row_name,input$column_name,input$by_column),{
    req(input$y_resp)

    row_name=input$row_name %||% TRUE
    column_name=input$column_name %||% TRUE
    by_column=input$by_column %||% FALSE

    file=read.csv(input$y_resp$datapath,row.names=if(row_name){1}else{NULL},header=column_name)
    if(by_column){
      file=as.data.frame(t(file))
    }
    names(file)=str_replace_all(names(file),' ','.')
    return(file)
    })
  data=eventReactive(c(input$y_resp,input$choosed_column,input$ref_multnom),{
    reac_vals$Series_n=input$column_index
    if(length(dim(raw_data()))>1){
      if(dim(raw_data())[2]>1){
        file=raw_data()[,names(raw_data()) %in% input$column_index]
        if(input$kernel=='Multinomial'){
          ref_col=file[,(names(file) %in% input$ref_multnom)]
          file=file[,!(names(file) %in% input$ref_multnom)]
          file=cbind(file,ref_col)
          names(file)[length(names(file))]=input$ref_multnom
        }
        if(is.null(dim(file))){
          file=data.frame(file)
          names(file)=input$column_index
          if(input$kernel=='Normal'){
            file=cbind(file,rep(0,length(file)))
            names(file)=c(input$column_index,'Variance')
            reac_vals$Series_n=c(reac_vals$Series_n,'Variance')
          }
        }
      }else{
        file=raw_data()[,1]
      }
    }
    return(file)
  })

  observeEvent(input$choosed_column,{
    removeModal()
    reac_vals$data_selected=TRUE
    })

  model=eventReactive(input$fit_model,{
    values_name=reac_vals$Series_n[!(reac_vals$Series_n %in% input$ref_multnom)]
    Series_n=dim(data())[2]-ifelse(input$kernel=='Multinomial',1,0)
    Time_length=dim(data())[1]
      aux_func=function(index_i){
        type_i=input[['type_'%>%paste0(index_i)]]
        order_i=input[['order_'%>%paste0(index_i)]]
        value_i=matrix(NA,Series_n,Time_length)

        for(index in 1:length(values_name)){
          label=if(input[['check_shared_F_'%>%paste0(index_i)]]){paste0(index_i)}else{paste0(index_i,'_',values_name[index])}
          pre_value_i=input[['value_'%>%paste0(label)]]
          if(is.na(pre_value_i)){
            row_name=input[['row_name_'%>%paste0(label)]] %||% TRUE
            column_name=input[['column_name_'%>%paste0(label)]] %||% TRUE
            by_column=input[['by_column_'%>%paste0(label)]] %||% FALSE

            file=read.csv(input[['x_resp_'%>%paste0(label)]]$datapath,row.names=if(row_name){1}else{NULL},header=column_name)
            if(by_column){
              file=as.data.frame(t(file))
            }
            names(file)=str_replace_all(names(file),' ','.')
            pre_value_i=file[[input[['custom_value_index_'%>%paste0(label)]]]]
          }
          value_i[index,]=pre_value_i
        }

        if(input[['check_offset_'%>%paste0(index_i)]]){
          if(input$kernel=='Multinomial'){
            label=paste0(index_i,'_',input$ref_multnom)
            pre_value_i=input[['value_'%>%paste0(label)]]
            if(is.na(pre_value_i)){
              row_name=input[['row_name_'%>%paste0(label)]] %||% TRUE
              column_name=input[['column_name_'%>%paste0(label)]] %||% TRUE
              by_column=input[['by_column_'%>%paste0(label)]] %||% FALSE

              file=read.csv(input[['x_resp_'%>%paste0(label)]]$datapath,row.names=if(row_name){1}else{NULL},header=column_name)
              if(by_column){
                file=as.data.frame(t(file))
              }
              names(file)=str_replace_all(names(file),' ','.')
              pre_value_i=file[[input[['custom_value_index_'%>%paste0(label)]]]]
            }
            for(index in c(1:Series_n)){
              value_i[index,]=value_i[index,]/pre_value_i
            }
          }
          value_i=log(value_i)
        }

        name_i=input[['name_'%>%paste0(index_i)]]
        D_i=input[['D_'%>%paste0(index_i)]]
        W_i=input[['W_'%>%paste0(index_i)]]
        m0_i=input[['m0_'%>%paste0(index_i)]]
        C0_i=input[['C0_'%>%paste0(index_i)]]
        create_block=if(type_i=='Polynomial'){gera_bloco_poly}else{gera_bloco_sazo}

        if(input[['check_shared_lat_'%>%paste0(index_i)]]){
          return(create_block(order_i,value=value_i,name=name_i,D=1/D_i,m0=m0_i,C0=C0_i,W=W_i))
        }else{
          aux_bloc=function(index_s){
            partial_value=value_i
            partial_value[-index_s,]=0
            create_block(order_i,value=partial_value,name=name_i %>% paste0('_',values_name[index_s]),D=1/D_i,m0=m0_i,C0=C0_i,W=W_i)
          }
          return(do.call(concat_bloco,lapply(1:length(values_name),aux_bloc)))
        }
      }
      structure=do.call(concat_bloco,lapply(reac_vals$par_index,aux_func))

      ajusta_modelo(data_out=data() %>% as.matrix,struture=structure,kernel=input$kernel)
    })

  observeEvent(input$fit_model,{
    reac_vals$current_fit_input=input %>% as.list
    reac_vals$current_par_index=reac_vals$par_index
    output$forecast_tab=renderUI({
      isolate({
        aux_func=function(index){
          if(input[['check_shared_F_'%>%paste0(index)]]){
            values_name=''
            show_values_name=''
            name_values_name=''
          }else{
            values_name=reac_vals$Series_n[!(reac_vals$Series_n %in% input$ref_multnom)]
            show_values_name=paste0('_',values_name)
            name_values_name=paste0(' (',values_name,')')
          }

          tabPanel(title=input[['name_'%>%paste0(index)]],
                   do.call(fluidRow,
                           lapply(X=1:length(values_name),
                                  FUN=function(index_s){
                                    num_input=numericInput('pred_value_'%>%paste0(index,show_values_name[index_s]),
                                                           paste0('Value',name_values_name[index_s],':'),
                                                           value=ifelse(length(values_name)==1,1,0))
                                    isolate({
                                      cond=!(show_values_name[index_s] %in% input$ref_multnom) | input[['check_offset_'%>%paste0(index)]]
                                    })
                                    toggleState(id='value_'%>%paste0(index,show_values_name[index_s]),
                                                condition=cond)
                                    num_input
                                  }
                           )
                   )
          )
        }

        column(8,
               panel('Aditional data for forecasting:',
                     do.call(
                       tabsetPanel,
                       lapply(reac_vals$par_index,aux_func)
                     )
               ),
               panel(numericInput('steps_ahead','Steps ahead:',min=1,value=1)),
               panel(plotlyOutput('forecast_plot')))
      })
    })
  })

  output$forecast_plot=renderPlotly({
      FF=array(0,c(0,dim(model()$FF)[2],input$steps_ahead))
      current_fit_input=reac_vals$current_fit_input
      values_name=current_fit_input$Series_n[!(current_fit_input$Series_n %in% current_fit_input$ref_multnom)]
      Series_n=dim(data())[2]-ifelse(current_fit_input$kernel=='Multinomial',1,0)

      for(index_i in reac_vals$current_par_index){
        order_i=if(current_fit_input[['type_'%>%paste0(index_i)]]=='Polynomial'){current_fit_input[['order_'%>%paste0(index_i)]]}else{2}
        value_i=matrix(NA,Series_n,input$steps_ahead)

        for(index in 1:length(values_name)){
          label=if(current_fit_input[['check_shared_F_'%>%paste0(index_i)]]){paste0(index_i)}else{paste0(index_i,'_',values_name[index])}
          pre_value_i=input[['pred_value_'%>%paste0(label)]]
          value_i[index,]=pre_value_i
        }
        if(current_fit_input[['check_offset_'%>%paste0(index_i)]]){
          if(current_fit_input$kernel=='Multinomial'){
            label=paste0(index_i,'_',input$ref_multnom)
            pre_value_i=current_fit_input[['pred_value_'%>%paste0(label)]]
            for(index in c(1:Series_n)){
              value_i[index,]=value_i[index,]/pre_value_i
            }
          }
          value_i=log(value_i)
        }

        if(current_fit_input[['check_shared_lat_'%>%paste0(index_i)]]){
          pre_FF=array(0,c(order_i,Series_n,input$steps_ahead))
          pre_FF[1,,]=value_i
        }else{
          aux_bloc=function(index_s){
            partial_value=value_i
            partial_value[-index_s,]=0
            pre_FF=array(0,c(order_i,Series_n,input$steps_ahead))
            pre_FF[1,,]=value_i
            return(pre_FF)
          }
          pre_FF=do.call(abind_aux_array,lapply(1:length(values_name),aux_bloc))
        }
        FF=abind_aux_array(FF,pre_FF)
      }
      predict(model(),t=input$steps_ahead,FF=FF,labels=names(data()))$plot
    })

  output$main_panel=renderUI({
    mainPanel(
      tabsetPanel(
        tabPanel('Data visualization',
                 plotlyOutput('prediction_plot')
                 ),
        tabPanel('Parameter selection',
               'Insert grid'),
        tabPanel('Forecasting',
                 uiOutput('forecast_tab')),
        tabPanel('Downloads',
                 'Insert download'),
        tabPanel('Help',
                 'Insert help')
      )
      )
  })

  reac_vals=reactiveValues(par_index=c(),data_selected=NULL,current_fit_input=list(),current_par_index=c(),Series_n=c())

  output$y_resp_input=renderUI({isolate(fileInput('y_resp', 'Select target data:'))})
  output$extra_param=renderUI({
    row=NULL
    if(input$kernel %||% 'None'=='Multinomial'){
      row=selectInput('ref_multnom','Select the column to use as base:',
                      choices=input$column_index)
    }
    return(row)
    })
  output$add_variable_input=renderUI({isolate(disabled(actionButton('add_variable','',icon = icon('plus'))))})
  output$fit_model=renderUI({disabled(actionButton('fit_model','Fit model'))})

  output$side_panel=renderUI({
      sidebarPanel(
        fluidRow("Model structure interface"),
        fluidRow(uiOutput('y_resp_input')),
        fluidRow(uiOutput('extra_param')),
        fluidRow(tabsetPanel(id='Var_tabs')),
        fluidRow(uiOutput('add_variable_input'),uiOutput('fit_model'))
      )
    })

  observeEvent(input$y_resp,{
    enable('add_variable')
  })
  observeEvent(c(input$add_variable,reac_vals$data_selected),{
    req(reac_vals$data_selected)
    index=as.character(input$add_variable+1)
    reac_vals$par_index=c(reac_vals$par_index,input$add_variable+1)
    create_type_ui(index,input,output,reac_vals,model)
    create_order_ui(index,input,output,reac_vals,model)
    create_name_ui(index,input,output,reac_vals,model)
    create_value_ui(index,input,output,reac_vals,model)
    create_m0_ui(index,input,output,reac_vals,model)
    create_C0_ui(index,input,output,reac_vals,model)
    create_D_ui(index,input,output,reac_vals,model)
    create_W_ui(index,input,output,reac_vals,model)
    create_custom_value_ui(index,input,output,reac_vals,model)
    create_delete_ui(index,input,output,reac_vals,model)
    create_show_value_ui(index,input,output,reac_vals,model)
    create_offset_button(index,input,output,reac_vals,model)
    create_full_ui(index,input,output,reac_vals,model)
    appendTab(
      inputId = 'Var_tabs',
      tab = tabPanel(title=textOutput('label_'%>%paste0(index)),
                     uiOutput('full_UI_'%>%paste0(index)),
                     value = 'tab_'%>%paste0(index)),
      select = TRUE
    )
    if(length(reac_vals$par_index)==1){
      enable('fit_model')
    }

  })

  output$prediction_plot=renderPlotly({
    req(reac_vals$data_selected)
    ref_data=data()
    if(input$kernel=='Normal'){
      ref_data=ref_data[,-2]
    }
    #req(input$column_index)
    if(as.numeric(input$fit_model)==0){
    t_last=ifelse(is.null(dim(ref_data)[1]),length(ref_data),dim(ref_data)[1])

    max_value=calcula_max(ref_data-min(ref_data))[[3]]+min(ref_data)
    min_value=-calcula_max(-(ref_data-max(ref_data)))[[3]]+max(ref_data)

    date=row.names(raw_data())
    time=c(1:t_last)

    pre_data=cbind(data.frame(time,date),ref_data)

    plot_data=(pre_data %>% pivot_longer(3:dim(pre_data)[2]))

    plt=ggplot(plot_data)+
      geom_point(aes(x=time,y=value,color=name,fill=name,date=date))+
      scale_fill_hue('',na.value=NA)+
      scale_color_hue('',na.value=NA)+
      scale_y_continuous(name='$y_t$')+
      scale_x_continuous('Time',labels=row.names(raw_data())[round(c(0:10)/(10/t_last))[-1]],breaks=round(c(0:10)/(10/t_last))[-1])+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90))+
      coord_cartesian(ylim=c(min_value,max_value))

    return(ggplotly(plt))
    }else{
      return(show_fit(model(),smooth = F,t_offset=1,dinamic=TRUE,labels=names(ref_data))$plot)
    }
  })
}

shinyApp(ui = final_ui, server = server)
