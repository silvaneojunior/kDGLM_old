create_order_ui=function(index,input,output,reac_vals,model){
  output[['block_order_'%>%paste0(index)]]=renderUI({
    numericInput(paste0('order_',index),
                 ifelse(input[['type_'%>%paste0(index)]]=='Polynomial','Order:','Period:'),
                 value=1,
                 min=1,
                 step=1
    )
  })}

create_type_ui=function(index,input,output,reac_vals,model){
  output[['block_type_'%>%paste0(index)]]=renderUI({
    selectInput('type_'%>%paste0(index),'Type:',choices=c('Polynomial','Harmonic'))
  })}

create_value_ui=function(index,input,output,reac_vals,model){
  output[['block_value_'%>%paste0(index)]]=renderUI({
    numericInput('value_'%>%paste0(index),'Value:',value=1)
  })}

create_name_ui=function(index,input,output,reac_vals,model){
  output[['block_name_'%>%paste0(index)]]=renderUI({
    textInput('name_'%>%paste0(index),'Name:',value=isolate(input[['name_'%>%paste0(index)]]) %||% 'Variable_'%>%paste0(index))
  })}

create_D_ui=function(index,input,output,reac_vals,model){
  output[['block_D_'%>%paste0(index)]]=renderUI({
    sliderInput('D_'%>%paste0(index),'Discount factor:',min=0,max=1,value=1)
  })}

create_W_ui=function(index,input,output,reac_vals,model){
  output[['block_W_'%>%paste0(index)]]=renderUI({
    numericInput('W_'%>%paste0(index),'W:',value=0,min=0)
  })}

create_m0_ui=function(index,input,output,reac_vals,model){
  output[['block_m0_'%>%paste0(index)]]=renderUI({
    numericInput('m0_'%>%paste0(index),'Prior mean:',value=0)
  })}

create_C0_ui=function(index,input,output,reac_vals,model){
  output[['block_C0_'%>%paste0(index)]]=renderUI({
    numericInput('C0_'%>%paste0(index),'Prior variance:',value=1)
  })}

create_custom_value_ui=function(index,input,output,reac_vals,model){
  output[['block_custom_value_'%>%paste0(index)]]=renderUI({
    actionButton('custom_value_'%>%paste0(index),'',icon=icon('cogs'))
  })
  output[['ui_custom_value_aux_col_'%>%paste0(index)]]=renderUI({
    req(input[['x_resp_'%>%paste0(index)]])

    row_name=input[['row_name_'%>%paste0(index)]] %||% TRUE
    column_name=input[['column_name_'%>%paste0(index)]] %||% TRUE
    by_column=input[['by_column_'%>%paste0(index)]] %||% FALSE

    file=read.csv(input[['x_resp_'%>%paste0(index)]]$datapath,row.names=if(row_name){1}else{NULL},header=column_name)
    if(by_column){
      file=as.data.frame(t(file))
    }
    selectInput('custom_value_index_'%>%paste0(index),'Select the column to insert:',choices=names(file))
  })
  observeEvent(input[['custom_value_'%>%paste0(index)]],{
    showModal(modalDialog(
      fluidRow(
        column(8,
               fileInput('x_resp_'%>%paste0(index), 'Select data:'),
               uiOutput('ui_custom_value_aux_col_'%>%paste0(index))),
        column(4,
               fluidRow(materialSwitch('by_column_'%>%paste0(index),'Variables on rows:')),
               fluidRow(materialSwitch('column_name_'%>%paste0(index),'Column names:',value=TRUE)),
               fluidRow(materialSwitch('row_name_'%>%paste0(index),'Row names:',value=TRUE)))),
      footer=fluidRow(actionButton('update_custom_value_'%>%paste0(index),'Select'),
                      modalButton('Cancel'))
      ))
    disable('update_custom_value_'%>%paste0(index))
  })
  observeEvent(input[['x_resp_'%>%paste0(index)]],{
    enable('update_custom_value_'%>%paste0(index))
    })
  observeEvent(input[['update_custom_value_'%>%paste0(index)]],{
    removeModal()
    row_name=input[['row_name_'%>%paste0(index)]] %||% TRUE
    column_name=input[['column_name_'%>%paste0(index)]] %||% TRUE
    by_column=input[['by_column_'%>%paste0(index)]] %||% FALSE

    file=read.csv(input[['x_resp_'%>%paste0(index)]]$datapath,row.names=if(row_name){1}else{NULL},header=column_name)
    if(by_column){
      file=as.data.frame(t(file))
    }
    updateNumericInput(inputId='value_'%>%paste0(index),value=file[[input[['custom_value_index_'%>%paste0(index)]]]])
  })
}

create_show_value_ui=function(index,input,output,reac_vals,model){
  output[['block_show_value_'%>%paste0(index)]]=renderUI({
    actionButton('show_value_'%>%paste0(index),'',icon=icon('area-chart'))
  })

  output[['plot_value_'%>%paste0(index)]]=renderPlotly({
    plot_lat_var(model(),input[['name_'%>%paste0(index)]],smooth=TRUE,cut_off=10,IC_prob=0.95,dinamic=TRUE,tranform_y=function(y){y})
    })

  observeEvent(input[['show_value_'%>%paste0(index)]],{
    showModal(modalDialog(
      plotlyOutput('plot_value_'%>%paste0(index))
    ))
  })
}

create_full_ui=function(index,input,output,reac_vals,model){
  output[['full_UI_'%>%paste0(index)]]=renderUI({
    column(12,fluidRow(column(4,uiOutput('block_name_' %>%paste0(index)),style = "width:300px;"),
                       column(4,uiOutput('block_check_offset_'%>%paste0(index)),style = "width:120px; padding-top:25px;"),
                     column(4 ,uiOutput('block_delete_'%>%paste0(index)),style = "width:50px; padding-top:25px;")),
         fluidRow(
           column(3,
                  fluidRow(uiOutput('block_type_'%>%paste0(index)),style = "width:120px;"),
                  fluidRow(uiOutput('block_order_'%>%paste0(index)),style = "width:120px; padding-top:20px;")),
           column(5,
                  fluidRow(uiOutput('block_D_'%>%paste0(index),style = "width:200px;")),
                  fluidRow(uiOutput('block_W_'%>%paste0(index),style = "width:200px;"))),
           column(2,
                  fluidRow(uiOutput('block_m0_'%>%paste0(index)),style = "width:120px;"),
                  fluidRow(uiOutput('block_C0_'%>%paste0(index)),style = "width:120px; padding-top:25px;"))),
         fluidRow(column(2,uiOutput('block_value_'%>%paste0(index)),style = "width:120px;"),
                  column(2,uiOutput('block_custom_value_'%>%paste0(index)),style = "width:120px; padding-top:25px;"),
                  column(2,uiOutput('block_show_value_'%>%paste0(index)),style = "width:120px; padding-top:25px;")),
         style='padding-right: 75px; padding-left: 50px;')
  })
  output[['label_'%>%paste0(index)]]=renderText(input[['name_'%>%paste0(index)]])
}

create_delete_ui=function(index,input,output,reac_vals,model){
  output[['block_delete_'%>%paste0(index)]]=renderUI({
    actionButton('delete_'%>%paste0(index),'',icon=icon('trash'))
  })
  observeEvent(input[['delete_'%>%paste0(index)]],{
    reac_vals$par_index=reac_vals$par_index[reac_vals$par_index!=index]
    output[['block_order_'%>%paste0(index)]]=NULL
    output[['block_type_'%>%paste0(index)]]=NULL
    output[['block_value_'%>%paste0(index)]]=NULL
    output[['block_name_'%>%paste0(index)]]=NULL
    output[['block_D_'%>%paste0(index)]]=NULL
    output[['block_W_'%>%paste0(index)]]=NULL
    output[['block_m0_'%>%paste0(index)]]=NULL
    output[['block_C0_'%>%paste0(index)]]=NULL
    output[['block_delete_'%>%paste0(index)]]=NULL
    output[['block_custom_value_'%>%paste0(index)]]=NULL
    output[['ui_custom_value_aux_col_'%>%paste0(index)]]=NULL
    output[['block_show_value_'%>%paste0(index)]]=NULL
    output[['plot_value_'%>%paste0(index)]]=NULL
    output[['full_UI_'%>%paste0(index)]]=NULL

    removeTab(
      inputId = 'Var_tabs',
      target = 'tab_'%>%paste0(index)
    )
    if(length(reac_vals$par_index)==0){
      disable('fit_model')
    }
    })
}

create_offset_button=function(index,input,output,reac_vals,model){
  output[['block_check_offset_'%>%paste0(index)]]=renderUI({
    checkboxInput('check_offset_'%>%paste0(index),'Is offset?')
  })
  observeEvent(input[['check_offset_'%>%paste0(index)]],{
    if(input[['check_offset_'%>%paste0(index)]]){
      updateSliderInput(inputId='D_'%>%paste0(index),value=1)
      updateNumericInput(inputId='W_'%>%paste0(index),value=0)
      updateNumericInput(inputId='m0_'%>%paste0(index),value=1)
      updateNumericInput(inputId='C0_'%>%paste0(index),value=0)
      updateSelectInput(inputId='type_'%>%paste0(index),selected='Polynomial')
      updateNumericInput(inputId='order_'%>%paste0(index),value=1)
    }else{
      updateNumericInput(inputId='m0_'%>%paste0(index),value=0)
      updateNumericInput(inputId='C0_'%>%paste0(index),value=1)
    }

    toggleState('D_'%>%paste0(index),!input[['check_offset_'%>%paste0(index)]])
    toggleState('W_'%>%paste0(index),!input[['check_offset_'%>%paste0(index)]])
    toggleState('m0_'%>%paste0(index),!input[['check_offset_'%>%paste0(index)]])
    toggleState('C0_'%>%paste0(index),!input[['check_offset_'%>%paste0(index)]])
    toggleState('type_'%>%paste0(index),!input[['check_offset_'%>%paste0(index)]])
    toggleState('order_'%>%paste0(index),!input[['check_offset_'%>%paste0(index)]])
    toggleState('show_value_'%>%paste0(index),!input[['check_offset_'%>%paste0(index)]])
  })
}

