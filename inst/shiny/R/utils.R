bold = function(x) strong(x, .noWS = "outside")
ital = function(x) em(x, .noWS = "outside")
link = function(s, href = s) a(s, href = href, .noWS = "outside")

mylink = function(text, href, .noWS = "outside", ...) {
  if(missing(href))
    href = text
  shiny::a(text, href = href, .noWS = .noWS, target = "_blank", ...)
}

errModal = function(..., html = FALSE) {
  args = list(...)
  if(length(args) == 1 && inherits(args[[1]], "condition"))
    mess = conditionMessage(args[[1]])
  else
    mess = paste(lapply(args, toString), collapse = "")
  if(html)
    mess = HTML(mess)

  showModal(modalDialog(mess, easyClose = TRUE))
}

parsePlotError = function(e) {
  msg = conditionMessage(e)
  if(grepl("reduce cex", msg))
    msg = "Too big for plot window!\n(You may still do simulations.)"
  msg
}


loadPed = function(file) {
  if(is.null(file))
    return()
  
  if(!file.exists(file)) 
    stop("File not found")
  
  df = read.table(file, header = TRUE, sep = "\t", colClasses = "character",
                  check.names = FALSE)
  names(df) = nms = tolower(names(df))

  cls = c("id", "fid", "mid", "sex")
  if(!all(cls %in% nms))
    stop("Column not found: ", toString(setdiff(cls, nms)))
  
  as.ped(df[cls])
}

checkSimInput = function(ped, ids, analysis, N) {
  if(is.null(ped)) 
    return("No pedigree indicated")
  if(length(ids) == 0) 
    return("No pedigree members indicated")
  if(!all(ids %in% labels(ped)))
    return(paste("Unknown ID label:", toString(setdiff(ids, labels(ped)))))
  if(analysis == "Sharing" && length(ids) == 1)
    return(paste("Sharing analysis is indicated, but only one individual:", toString(ids)))
  if(analysis == "Autozygosity" && any((inbr <- inbreeding(ped, ids)) == 0))
    return(paste("Autozygosity analysis indicated, but some individuals are not inbred: ", toString(names(inbr)[inbr == 0])))
  if(!ibdsim2:::isCount(N))
    return("Number of simulations must be a positive integer")
  if(N < 5)
    return("Number of simulations is too low (min 5)")
  if(N > 5000)
    return("Number of simulations is too high (max 5000)")
  
  "ok"
}

myMainPanel = function(...) {
  div(class = "col-sm-6 col-lg-8", role = "main", ...)
}

mySidebarPanel = function(...) {
  div(class = "col-sm-3 col-lg-2", 
      tags$form(class = "well", role = "complementary", ...))
}

getMapLength = function(map, unit, chrom) {
    if(unit == "mb") physRange(map) 
    else if(chrom == "X") mapLen(map, sex = "female")
    else mean(mapLen(map))
}