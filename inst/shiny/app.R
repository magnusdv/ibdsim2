suppressMessages(suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinyWidgets)
  library(pedtools)
  library(ribd)
  library(ibdsim2)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(glue)
  library(zip)
}))


.VERSION = packageDescription("ibdsim2")$Version

# Variables used in both side panels
.MODNAMES = list(HTML("&#x03C7;<sup>2</sup>"), "Haldane")
.MODVALS  = c("chi", "haldane")


ui = fluidPage(
  includeCSS("www/custom.css"),
  tags$head(includeHTML(system.file("shiny/www/GA.html", package = "ibdsim2"))),
  tags$head(tags$script(src = "scripts.js")),
  
  useShinyjs(),
  useBusyIndicators(),
  
  tags$head(tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Lobster&display=swap")),
  
  # tags$div(id = "banner",
  #       p(id="big-text", "Major app update!"),
  #       p("Check out the ", mylink("NEWS", href="https://github.com/magnusdv/ibdsim2/blob/master/NEWS.md", 
  #                                 style = "font-weight:bold;")),
  #       #p(id="small-text", "The old version still available ", 
  #       #  mylink("here", href="https://magnusdv.shinyapps.io/ibdsim2-14/"))
  # ),

  # Application title
  h2(id = "title-h2", "IBD sharing by family members"),
  
  p(style = "margin-bottom: 4px", bold("Purpose: "),
"Estimate and visualise distributions of genomic segments shared identical-by-descent (IBD) between related individuals, 
or within inbred individuals (autozygosity). Recombination is simulated down through the pedigree, using detailed, sex-specific crossover rates for the human genome (",
mylink("Halldorsson et al., 2019", "https://doi.org/10.1126/science.aau1043"), ")."),

  p(style = "margin-bottom: 4px", bold("More information: "),
    "This program is a frontend for the R package ", mylink("ibdsim2", "https://github.com/magnusdv/ibdsim2"), 
    ", which is part of the ", mylink("pedsuite", "https://magnusdv.github.io/pedsuite"), " ecosystem for pedigree analysis.", 
    "Details about the simulations can be found in the package documentation, and in the book ",
    mylink("Pedigree analysis in R", "https://www.elsevier.com/books/pedigree-analysis-in-r/vigeland/978-0-12-824430-2"), 
    ". Please cite this book if you use the app in your work."), 

# Widgets --------------------------------------------------------------
fluidRow(
  
  # Left sidebar
  mySidebarPanel( # style = "padding-top: 5px; padding-bottom:5px",
    h4("Pedigree 1"),
    selectizeInput("builtin1", "Built-in pedigree", selected = "Half-sibs, maternal", 
                   choices = c(Choose = "", names(BUILTIN_PEDS)), size = 10),
    fileInput("loadped1", "Load ped file", buttonLabel = icon("folder-open"),
              accept = c(".ped", ".txt"), width = "100%", placeholder = NULL),
   textInput("ids1", "Individuals", value = "", width = "100%"),
   textInput("label1", "Label", value = "Ped 1", width = "100%"),

    hr(style = "border-top: 1px solid #000000; margin-top: 5px; margin-bottom: 0px"),
    
    div(id = "opt1", h3("PARAMETERS", id = "paramHeading"),
      radioButtons("chrom1", "Chromosome", choices = c("1 - 22" = "aut", "X" = "X"), 
                   selected = "aut", inline = TRUE),
      radioButtons("model1", "Crossover model", choiceNames = .MODNAMES, choiceValues = .MODVALS, 
                   selected = "chi", inline = TRUE),
      radioButtons("sexspec1", "Sex-specific map", choices = c("On", "Off"), 
                   selected = "On", inline = TRUE),
      fluidRow(
        column(6, class = "leftinput", 
               numericInput("cutoff1", "Cutoff", value = 0, min = 0, step = 1)),
        column(6, class = "rightinput", 
               numericInput("seed1", "Seed", value = 123, min = 1, step = 1)),
      ),
    ),
    
    hr(style = "border-top: 1px solid #000000; margin-top: 10px; margin-bottom: 20px"),
    
    fluidRow(
      column(9, actionButton("simulate1", "Simulate!", width = "100%", class = "btn btn-primary btn-lg")),
      column(3, style = "padding-left:5px;", uiOutput("icon1"))
    ),
   # actionButton("rcode1", "R code"),
  ),
  
  # Middle region: Plots
  myMainPanel(
    fluidRow(
      column(6, align = "center", plotOutput("pedplot1", height = "295px")),
      column(6, align = "center", plotOutput("pedplot2", height = "295px"))
    ),
    plotOutput("ibdplot", width = "100%"),
    absolutePanel(id = "origin0panel",
      bottom = 0, left = 0, width = "auto", draggable = FALSE, style = "z-index:1000",
      prettySwitch("orig0", label = "Origin", value = FALSE, slim = TRUE, status = "default")
    )
  ),
  
  # Right sidebar
  mySidebarPanel( # style = "padding-top: 5px; padding-bottom:5px",
    h4("Pedigree 2"),
    selectizeInput("builtin2", "Built-in pedigree", selected = "", 
                   choices = c(Choose = "", names(BUILTIN_PEDS)), size = 10),
    fileInput("loadped2", "Load ped file", buttonLabel = icon("folder-open"),
              accept = c(".ped", ".txt"), width = "100%", placeholder = NULL),
    textInput("ids2", "Individuals", value = "", width = "100%"),
    textInput("label2", "Label", value = "Ped 2", width = "100%"),
    
    hr(style = "border-top: 1px solid #000000; margin-top: 5px; margin-bottom: 0px"),
    
    div(id = "opt2", h3("PARAMETERS", id = "paramHeading"),
      radioButtons("chrom2", "Chromosome", choices = c("1 - 22" = "aut", "X" = "X"), 
                   selected = "aut", inline = TRUE),
      radioButtons("model2", "Crossover model", choiceNames = .MODNAMES, choiceValues = .MODVALS, 
                   selected = "chi", inline = TRUE),
      radioButtons("sexspec2", "Sex-specific map", choices = c("On", "Off"), 
                   selected = "On", inline = TRUE),
      fluidRow(
        column(6, class = "leftinput", 
               numericInput("cutoff2", "Cutoff", value = 0, min = 0, step = 1)),
        column(6, class = "rightinput", 
               numericInput("seed2", "Seed", value = 123, min = 1, step = 1)),
      ),
    ),
    
    hr(style = "border-top: 1px solid #000000; margin-top: 10px; margin-bottom: 20px"),
    fluidRow(
      column(9, actionButton("simulate2", "Simulate!", width = "100%", class = "btn btn-primary btn-lg")),
      column(3, style = "padding-left:5px;", uiOutput("icon2"))
    ),
  ),
),
  
# Bottom panel
fluidRow(
  column(6, wellPanel(id = "bottomwell1", style = "position:relative",
    div(style = "position:absolute; right: 10px; top: 15px; z-index: 1000",
        downloadButton("download", "Download data", class="btn btn-success btn-sm")),
    fluidRow(
      column(4, 
        h4("Settings"),
        numericInput("nsims", "Number of sims", value = 500, min = 5, max = 5000),
      ),
      column(8, 
        radioButtons("unit", "Length unit", selected = "cm", inline = TRUE, 
                       choices = c("cM" = "cm", "Mb" = "mb"), width = "100%"),
        radioButtons("analysis", "Analysis", selected = "Sharing", inline = TRUE,
                     choices = c("Sharing", "Autozygosity"), width = "100%"),
      ),
  ))),
  column(6, wellPanel(id = "bottomwell2",
    fluidRow(
      column(8, h4("Observed data"),
        fluidRow(
          column(6, numericInput("obs-total", "Total length", value = "")),
          column(6, numericInput("obs-nseg", "Count", value = "")),
      )),
      column(4, textAreaInput("obs-segs", "Segments", value = "", rows = 4)),
    ))),
  ),

  p(style = "font-size:small", "This is version", .VERSION, "of ibdsim2 (",
  mylink("changelog", "https://github.com/magnusdv/ibdsim2/blob/master/NEWS.md"), ").",
  "For bug reports, feature requests, or other comments, please file an issue at ", 
  mylink("https://github.com/magnusdv/ibdsim2/issues"), "."),
)


# Server logic
server = function(input, output, session) {

  observeEvent(input$browserClosed, stopApp())
  
  ped1 = reactiveVal(NULL)
  ped2 = reactiveVal(NULL)
  
  ids1 = reactive(setdiff(trimws(strsplit(input$ids1, ",")[[1]]), ""))
  ids2 = reactive(setdiff(trimws(strsplit(input$ids2, ",")[[1]]), ""))
  
  observeEvent(input$builtin1, {
    choice = req(input$builtin1)
    ped1(BUILTIN_PEDS[[choice]])
    updateTextInput(session, "ids1", value = toString(DEFAULT_IDS[[choice]]))
  })
  
  observeEvent(input$builtin2, {
    choice = req(input$builtin2)
    ped2(BUILTIN_PEDS[[choice]])
    updateTextInput(session, "ids2", value = toString(DEFAULT_IDS[[choice]]))
  })

  observeEvent(input$loadped1, {
    x = tryCatch(loadPed(input$loadped1$datapath), 
                   error = errModal, warning = errModal)
    ped1(req(x))
    updateTextInput(session, "ids1", value = "")
    isolate(updateSelectizeInput(session, "builtin1", selected = ""))
  })
  
  observeEvent(input$loadped2, {
    ped = tryCatch(loadPed(input$loadped2$datapath),
                   error = errModal, warning = errModal)
    ped2(req(ped))
    updateTextInput(session, "ids2", value = "")
    isolate(updateSelectizeInput(session, "builtin2", selected = ""))
  })

  observeEvent(input$chrom1, {
    if(input$chrom1 == "X") {
      updateRadioButtons(session, "sexspec1", selected = "On")
      disable("sexspec1")
    }
    else
      enable("sexspec1")
  })
  
  observeEvent(input$chrom2, {
    if(input$chrom2 == "X") {
      updateRadioButtons(session, "sexspec2", selected = "On")
      disable("sexspec2")
    }
    else
      enable("sexspec2")
  })
  
  map1 = reactive({
    chr = switch(input$chrom1, aut = 1:22, X = 23)
    unif = tolower(input$unit) == "cm"
    sexspec = if(input$chrom1 == "X") TRUE else input$sexspec1 == "On"
    loadMap("decode19", chrom = chr, uniform = unif, sexAverage = !sexspec)
  })

  map2 = reactive({
    chr = switch(input$chrom2, aut = 1:22, X = 23)
    unif = tolower(input$unit) == "cm"
    sexspec = if(input$chrom2 == "X") TRUE else input$sexspec2 == "On"
    loadMap("decode19", chrom = chr, uniform = unif, sexAverage = !sexspec)
  })
  
  maplen1 = reactive(getMapLength(map1(), input$unit, input$chrom1))
  maplen2 = reactive(getMapLength(map2(), input$unit, input$chrom2))
  
# Simulations -------------------------------------------------------------

  sim1 = reactiveVal(NULL)
  sim2 = reactiveVal(NULL)
 
  # Reset if anything changes
  observe({ped1(); ids1(); map1(); input$model1; input$nsims; input$seed1; input$analysis; sim1(NULL); enable("simulate1")})
  observe({ped2(); ids2(); map2(); input$model2; input$nsims; input$seed2; input$analysis; sim2(NULL); enable("simulate2")})
  
  # Icons
  output$icon1 = renderUI(icon(name = if(is.null(sim1())) "arrow-left" else "check"))
  output$icon2 = renderUI(icon(name = if(is.null(sim2())) "arrow-left" else "check"))
  
  # Simulate!
  observeEvent(input$simulate1, {
    chk = checkSimInput(ped1(), ids1(), input$analysis, input$nsims)
    if(chk != "ok")
      return(errModal(chk))
    disable("simulate1")
    sims = ibdsim(ped1(), N = input$nsims, ids = ids1(), map = map1(), 
                  model = input$model1, seed = input$seed1, verbose = FALSE) 
    sim1(sims)
  })
  
  observeEvent(input$simulate2, {
    chk = checkSimInput(ped2(), ids2(), input$analysis, input$nsims)
    if(chk != "ok")
      return(errModal(chk))
    disable("simulate2")
    sims = ibdsim(ped2(), N = input$nsims, ids = ids2(), map = map2(), 
                  model = input$model2, seed = input$seed2, verbose = FALSE)
    sim2(sims)
  })
  
  segmentData1 = reactive(getSegmentData(sim1(), analysis = input$analysis, cutoff = input$cutoff1, unit = input$unit))
  segmentData2 = reactive(getSegmentData(sim2(), analysis = input$analysis, cutoff = input$cutoff2, unit = input$unit))
  

# Observed data -----------------------------------------------------------

  observedSegs = reactive({
    lenStr = input$`obs-segs` |> strsplit("\n") |> unlist() |> strsplit(",") |> unlist() |> trimws()
    lenStr = lenStr[lenStr != ""]
    lens = suppressWarnings(as.numeric(lenStr))
    if(anyNA(lens))
      return(errModal("Non-numeric segment length: ", lenStr[is.na(lens)]))
    lens
  })
  
  observeEvent(input$`obs-segs`, {
    lens = observedSegs()
    if(!length(lens)) {
      enable("obs-nseg"); enable("obs-total")
      updateNumericInput(session, "obs-nseg", value = "")
      updateNumericInput(session, "obs-total", value = "")
    }
    else {
      updateNumericInput(session, "obs-nseg", value = length(lens))
      updateNumericInput(session, "obs-total", value = sum(lens))
      disable("obs-nseg"); disable("obs-total")
    }
  })
  
  observed = reactive({
    nseg = input$`obs-nseg`
    total = input$`obs-total`
    if(is.na(nseg) || is.na(total))
      return(NULL)
    list(nseg = nseg, total = total, mean = total/nseg, lengths = observedSegs())
  })
  
# Plots ----------------------------------------------------------
  
  # Red and blue used consistently for the two pedigrees
  COLS = c(2, 4)
  
  output$pedplot1 = renderPlot({
    ped = req(ped1())
    lab = input$label1
    plotWidth = session$clientData$output_pedplot1_width
    mar = calculateMargin(plotWidth)
    isBuiltin = input$builtin1 != ""
    
    tryCatch(
      plotped(ped, ids = ids1(), col = COLS[1], title = lab, margin = mar, 
              straightlegs = isBuiltin),
      error = function(e) {
        plot.new(); box(which = "outer", col = 1); title(lab); 
        text(x = 0.5, y = 0.6, parsePlotError(e), cex = 1.1, col = 2)
      })
  }, execOnResize = TRUE)

  output$pedplot2 = renderPlot({
    ped = req(ped2())
    lab = input$label2
    plotWidth = session$clientData$output_pedplot2_width
    mar = calculateMargin(plotWidth)
    isBuiltin = input$builtin2 != ""
    
    tryCatch(
      plotped(ped, ids = ids2(), col = COLS[2], title = lab, margin = mar,
              straightlegs = isBuiltin),
      error = function(e) {
        plot.new(); box(which = "outer", col = 1); title(lab); 
        text(x = 0.5, y = 0.6, parsePlotError(e), cex = 1.1, col = 2)
      })
  }, execOnResize = TRUE)
  
  output$ibdplot = renderPlot({
    labs = c(input$label1, input$label2)
    segData = list(segmentData1(), segmentData2())
    cols = COLS
    names(segData) = names(cols) = labs
    
    isnull = sapply(segData, is.null)
    skip = isnull | labs == ""
    
    for(i in 1:2) if(!isnull[i] && labs[i] == "")
      return(errModal("Please specify a label for pedigree ", i))
    
    if(!any(skip) && labs[1] == labs[2])
      return(errModal("Labels cannot be equal"))
    
    req(!isnull)  # return if both empty
    
    # Map length (for percentages)
    if(!any(skip)) maplen = if(maplen1() == maplen2()) maplen1() else NULL
    else maplen = if(skip[1]) maplen2() else maplen1()
    
    g = generateIbdPlot(segData[!skip], input$analysis, cols = cols[!skip],
                        unit = input$unit, observed = observed(), 
                        orig0 = input$orig0, maplen = maplen)
    suppressWarnings(print(g))
  })
  

# Generate R code----------------------------------------------------------

  # codeTxt = reactiveVal(NULL)
  # 
  # # Render in modal dialog, created with createCodeModal when pressing rcode button (see below)
  # output$showCode = renderText(req(codeTxt()))
  # 
  # output$saveCode = downloadHandler(
  #   filename = "ibdsim.R",
  #   content = function(con) {
  #     cat(codeTxt(), file = con)
  #     removeModal()
  #   },
  #   contentType = "text/plain"
  # )
  # 
  # observeEvent(input$rcode1, {
  #   code = generateCode(ped = ped1(),
  #                       ids = ids1(),
  #                       chrom = input$chrom1, model = input$model1, 
  #                       sexspec = input$sexspec1, cutoff = input$cutoff1, 
  #                       analysis = input$analysis, 
  #                       unit = input$unit, nsims = input$nsims, seed = input$seed1)
  #   codeTxt(code)
  #   showModal(createCodeModal())
  # })
  
# Download data -----------------------------------------------------------

  allParams1 = reactive(list(
    ped = ped1(), label = input$label1, builtin = input$builtin1, ids = ids1(), 
    loadped = input$loadped1$name, chrom = input$chrom1, model = input$model1, 
    sexspec = input$sexspec1, cutoff = input$cutoff1, analysis = input$analysis, 
    unit = input$unit, nsims = input$nsims, seed = input$seed1)
  )
  
  allParams2 = reactive(list(
    ped = ped2(), label = input$label2, builtin = input$builtin2, ids = ids2(),
    loadped = input$loadped2$name, chrom = input$chrom2, model = input$model2,
    sexspec = input$sexspec2, cutoff = input$cutoff2, analysis = input$analysis, 
    unit = input$unit, nsims = input$nsims, seed = input$seed2)
  )
  
  output$download = downloadHandler(
    filename = "ibdsim2-output.zip",
    content = function(con) {
      tmpdir = tempdir()
      files = saveData(segmentData1(), segmentData2(), params1 = allParams1(), 
                       params2 = allParams2(), version = .VERSION, tmpdir = tmpdir)
    
      if(!length(files)) return(errModal("No data to save"))
      zip::zip(con, files, root = tmpdir)
    }, 
    contentType = "application/zip"
  )
  
}

# Run the application 
suppressMessages(suppressPackageStartupMessages(
  shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
))
