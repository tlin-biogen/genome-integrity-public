# Purpose: Shiny app for ddPCR analysis of genome integrity
library(shiny)
library(tidyverse)
library(purrr)
library(polynom)
library(reactlog)
library(DT) # will mask dataTableOutput and renderDataTable

source("ddPCR_genome_integrity_global.R")

ui <- fluidPage(
  titlePanel(
    "2D-ddPCR data analysis tool"
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("datafile", "Please upload a CSV file", accept = c(".csv", ".tsv")),
      downloadButton("download", label = "Download Results", class = "btn-block"),
      width = 3
    ),
    mainPanel(
      # output 1
      titlePanel("1. Raw data"),
      DT::dataTableOutput("preview1"),
      titlePanel("2. Data with estimate"),
      DT::dataTableOutput("preview2"),
      width = 8
    )
  )
)


server <- function(input, output, session) {

  # Upload data
  raw <- reactive({

    req(input$datafile)

    infile <- read.csv(input$datafile$datapath, row.names = NULL, check.names = FALSE) 
    infile <- process_csv(infile)
    
    colnames(infile) <- toupper(colnames(infile))

    infile <- infile %>% select(WELL, SAMPLE, CONCENTRATION, TARGET, POSITIVES, NEGATIVES, starts_with("CH1"), ACCEPTEDDROPLETS)
    infile
  })

  output$preview1 <- DT::renderDataTable(raw())

  # Clean and estimate
  cleaned <- reactive({

    req(raw())

    out <- raw()

    out <- out %>% mutate(
      lambda = -log(`CH1-CH2-` / ACCEPTEDDROPLETS),

      m = floor(ACCEPTEDDROPLETS * lambda)
    )

    # Create a new id to allow for merging
    out$id <- seq(1, nrow(out), 1)

    subout <- out %>% filter(`CH1+CH2-` != 0 & `CH1-CH2+` != 0 & `CH1-CH2-` != 0)

    subout <- subout %>% mutate(
      prob1 = list(m, ACCEPTEDDROPLETS, `CH1+CH2-`, rep(20, length(m))) %>% purrr::pmap_dbl(.f = polysolver),
      prob2 = list(m, ACCEPTEDDROPLETS, `CH1-CH2+`, rep(20, length(m))) %>% purrr::pmap_dbl(.f = polysolver),
      
      est = (1 - prob1 - prob2) * 100
    )

    out <- left_join(out, subout %>% select(id, est), by = "id") %>%
      mutate(ESTIMATED_PERCENTAGE = case_when(
        toupper(CONCENTRATION) == "NO CALL" ~ "NO CALL",
        `CH1-CH2-` == 0 | `CH1+CH2-` == 0 | `CH1-CH2+` == 0 ~ "NA",
        TRUE ~ format(round(est, digits = 2), nsmall = 2)
      )) %>%
      group_by(SAMPLE) %>%
      arrange(SAMPLE, id) %>%
      mutate(
        SAME_EST_WITHIN_SAMPLE = case_when(
          is.na(est[1]) & is.na(est[2]) ~ NA,
          is.na(est[1]) != is.na(est[2]) ~ FALSE,
          TRUE ~ abs((est[1] - est[2])) < 1e-8
        )
      ) %>%
      ungroup() %>%
      select(WELL, SAMPLE, CONCENTRATION, TARGET, POSITIVES, NEGATIVES, starts_with("CH1"), ACCEPTEDDROPLETS, ESTIMATED_PERCENTAGE, SAME_EST_WITHIN_SAMPLE)

    out
  })

  output$preview2 <- DT::renderDataTable(cleaned() %>% DT::datatable() %>% 
                                           DT::formatStyle("SAME_EST_WITHIN_SAMPLE", backgroundColor = styleEqual(c(TRUE, FALSE), c("gray", "yellow"))))

  # Download
  output$download <- downloadHandler(
    filename = function() {
      "output.csv"
    },
    content = function(file) {
      write.csv(cleaned(), file)
    }
  )
}


shinyApp(ui = ui, server = server)