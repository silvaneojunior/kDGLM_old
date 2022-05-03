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

# Misc stuff
library(hms)
library(beepr)
library(latex2exp)

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

server=function(input, output) {

  modal_prompt=modalDialog(
    #numericInput('column_index', 'Select the index of the column to insert:', value = 1,min=1,max=,step=1),
    fluidRow(
      column(8,
             uiOutput('modal_content'),
      column(4,
             fluidRow(materialSwitch('by_column','Variables on rows:')),
             fluidRow(materialSwitch('column_name','Column names:',value=TRUE)),
             fluidRow(materialSwitch('row_name','Row names:',value=TRUE))))),
    footer = actionButton('choosed_column','Select')
  )

  output$modal_content=renderUI({
             selectInput('column_index','Select the column to insert:',choices=names(raw_data()))
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
    return(file)
    })
  data=eventReactive(c(input$y_resp,input$choosed_column),{
    if(length(dim(raw_data()))>1){
      if(dim(raw_data())[2]>1){
        file=raw_data()[[input$column_index]]
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
      blocks_list=reac_vals$par_index
      index_i=blocks_list[1]
      type_i=input[['type_'%>%paste0(index_i)]]
      order_i=input[['order_'%>%paste0(index_i)]]
      value_i=input[['value_'%>%paste0(index_i)]]
      name_i=input[['name_'%>%paste0(index_i)]]
      D_i=input[['D_'%>%paste0(index_i)]]
      W_i=input[['W_'%>%paste0(index_i)]]
      m0_i=input[['m0_'%>%paste0(index_i)]]
      C0_i=input[['C0_'%>%paste0(index_i)]]

      if(is.na(value_i)){
        row_name=input[['row_name_'%>%paste0(index_i)]] %||% TRUE
        column_name=input[['column_name_'%>%paste0(index_i)]] %||% TRUE
        by_column=input[['by_column_'%>%paste0(index_i)]] %||% FALSE

        file=read.csv(input[['x_resp_'%>%paste0(index_i)]]$datapath,row.names=if(row_name){1}else{NULL},header=column_name)
        if(by_column){
          file=as.data.frame(t(file))
        }
        value_i=file[[input[['custom_value_index_'%>%paste0(index_i)]]]]
      }

      create_block=if(type_i=='Polynomial'){gera_bloco_poly}else{gera_bloco_sazo}
      structure=create_block(order_i,value=value_i,name=name_i,D=1/D_i,m0=m0_i,C0=C0_i,W=W_i)
      for(index_i in blocks_list[-1]){
        type_i=input[['type_'%>%paste0(index_i)]]
        order_i=input[['order_'%>%paste0(index_i)]]
        value_i=input[['value_'%>%paste0(index_i)]]
        name_i=input[['name_'%>%paste0(index_i)]]
        D_i=input[['D_'%>%paste0(index_i)]]
        W_i=input[['W_'%>%paste0(index_i)]]
        m0_i=input[['m0_'%>%paste0(index_i)]]
        C0_i=input[['C0_'%>%paste0(index_i)]]

        if(is.na(value_i)){

          row_name=input[['row_name_'%>%paste0(index_i)]] %||% TRUE
          column_name=input[['column_name_'%>%paste0(index_i)]] %||% TRUE
          by_column=input[['by_column_'%>%paste0(index_i)]] %||% FALSE

          file=read.csv(input[['x_resp_'%>%paste0(index_i)]]$datapath,row.names=if(row_name){1}else{NULL},header=column_name)
          if(by_column){
            file=as.data.frame(t(file))
          }
          value_i=file[[input[['custom_value_index_'%>%paste0(index_i)]]]]
        }

        create_block=if(type_i=='Polynomial'){gera_bloco_poly}else{gera_bloco_sazo}

        structure=concat_bloco(structure,create_block(order_i,value=value_i,name=name_i,D=1/D_i,m0=m0_i,C0=C0_i,W=W_i))
      }

      ajusta_modelo(data_out=data(),struture=structure)
    })

  observeEvent(input$fit_model,{
    output$forecast_tab=renderUI({
      column(8,
        panel('Aditional data for forecasting:',
                 do.call(
                   tabsetPanel,
                   map(names(model()$names),function(var){
                     tabPanel(var,textInput('pred_value_'%>%paste0(var),label='Variable value:',value=1))
                     })
                   )
                 ),
        panel(numericInput('steps_ahead','Steps ahead:',min=1,value=1)),
        panel(plotlyOutput('forecast_plot')))
    })
  })
  output$forecast_plot=renderPlotly({
      FF=matrix(0,0,input$steps_ahead)
      for(var in names(model()$names)){
        value=input[['pred_value_'%>%paste0(var)]]
        value=value %>% str_remove_all(' ') %>% str_split(',')
        value=value[[1]] %>% as.numeric
        if(length(value)==1){
          value=rep(value,input$steps_ahead)
        }
        if(length(value)<input$steps_ahead){
          value=c(value,rep(1,input$steps_ahead-length(value)))
        }else{if(length(value)>input$steps_ahead){
          value=value[1:input$steps_ahead]
        }
        }

        n_inputs=length(model()$names[[var]])
        pre_FF=matrix(0,n_inputs,input$steps_ahead)
        pre_FF[1,]=value
        FF=rbind(FF,pre_FF)
      }
      predict(model(),t=input$steps_ahead,FF=FF)$plot
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

  reac_vals=reactiveValues(par_index=c(),data_selected=NULL)

  output$y_resp_input=renderUI({isolate(fileInput('y_resp', 'Select target data:'))})
  output$add_variable_input=renderUI({isolate(disabled(actionButton('add_variable','',icon = icon('plus'))))})
  output$fit_model=renderUI({disabled(actionButton('fit_model','Fit model'))})

  output$side_panel=renderUI({
      sidebarPanel(
        fluidRow("Model structure interface"),
        fluidRow(uiOutput('y_resp_input')),
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
    #req(input$column_index)
    if(as.numeric(input$fit_model)==0){
    t_last=length(data())

    fill_list=c('black')
    names(fill_list)=c('Observed values')
    color_list=c('black')
    names(color_list)=c('Observed values')

    max_value=calcula_max(data()-min(data()))[[3]]+min(data())
    min_value=-calcula_max(-(data()-max(data())))[[3]]+max(data())

    date=row.names(raw_data())

    plt=ggplot()+geom_point(aes(x=c(1:t_last),y=data(),color='Observed values',fill='Observed values',date=date))+
      scale_fill_manual('',na.value=NA,values=fill_list)+
      scale_color_manual('',na.value=NA,values=color_list)+
      scale_y_continuous(name='$y_t$')+
      scale_x_continuous('Time',labels=row.names(raw_data())[round(c(0:10)/(10/t_last))[-1]],breaks=round(c(0:10)/(10/t_last))[-1])+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90))+
      coord_cartesian(ylim=c(min_value,max_value))

    return(ggplotly(plt))
    }else{
      return(show_fit(model(),smooth = T,t_offset=1,dinamic=TRUE)$plot)
    }
  })
}

shinyApp(ui = final_ui, server = server)
