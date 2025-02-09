##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing LoadingBoard")

LoadingInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shinyBS::tipify( shiny::actionLink(ns("module_info"), "Tutorial", icon = shiny::icon("youtube")),
               "Show more information about this module.")
        ## shiny::uiOutput(ns("inputsUI"))
        ## shiny::uiOutput(ns("socialButtons"))
    )
}

LoadingUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Datasets",uiOutput(ns("pgxtable_UI"))),
            shiny::tabPanel("Upload data",uiOutput(ns("upload_UI")))
            ## shiny::tabPanel("Community forum",uiOutput(ns("forum_UI")))
        )
    )
}

LoadingBoard <- function(input, output, session, pgx_dir=PGX.DIR, 
                         limits = c("samples"=1000,"comparisons"=20,
                                    "genes"=20000, "genesets"=10000,
                                    "datasets"=10),
                         enable_upload = TRUE,
                         enable_delete = TRUE,
                         enable_save = TRUE,
                         enable_userdir = TRUE,                         
                         authentication="none")
{
    ns <- session$ns ## NAMESPACE
    dbg("[LoadingBoard] >>> initializing LoadingBoard...")

    loadedDataset <- shiny::reactiveVal(FALSE)
    
    SHOWSPLASH=TRUE
    ## SHOWSPLASH=FALSE

    message("[LoadingBoard] in.shinyproxy = ",in.shinyproxy())    
    message("[LoadingBoard] SHINYPROXY_USERNAME = ",Sys.getenv("SHINYPROXY_USERNAME"))
    message("[LoadingBoard] SHINYPROXY_USERGROUPS = ",Sys.getenv("SHINYPROXY_USERGROUPS"))
    message("[LoadingBoard] authentication = ",authentication)

    auth <- NULL   ## shared in module
    if(authentication == "password") {
        auth <- shiny::callModule(
            PasswordAuthenticationModule, "auth",
            credentials.file = "CREDENTIALS")
    } else if(authentication == "firebase") {
        auth <- shiny::callModule(
                           FirebaseAuthenticationModule, "auth")
    } else if(authentication == "register") {
        auth <- shiny::callModule(
            RegisterAuthenticationModule, "auth",
            register.file = "../logs/register.log")
        ##} else if(authentication == "shinyproxy" && in.shinyproxy()) {
    } else if(authentication == "shinyproxy") {        
        username <- Sys.getenv("SHINYPROXY_USERNAME")
        ##email <- Sys.getenv("SHINYPROXY_EMAIL")        
        auth <- shiny::callModule(NoAuthenticationModule, "auth",
                                  username=username, email=username)
    } else {
        ## none
        auth <- shiny::callModule(NoAuthenticationModule, "auth")
    } 

    dbg("[LoadingBoard] names.auth = ",names(auth))
    
    ##-----------------------------------------------------------------------------
    ## Description
    ##-----------------------------------------------------------------------------
    description = "<b>Omics Playground</b> is a self-service bioinformatics platform for interactive analysis, visualization and interpretation of transcriptomics and proteomics data. Life scientists can easily perform complex data analysis and visualization without coding, and significantly reduce the time to discovery."
    
    output$description <- shiny::renderUI(shiny::HTML(description))

    shiny::observeEvent( input$module_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Data View Board</strong>"),
            shiny::HTML(module_infotext),
            easyClose = TRUE, size="l" ))
    })

    module_infotext =paste0(
        'The platform starts running from the <strong>Home panel</strong>. This panel shows the available datasets within the platform. The table reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date.

<br><br><b>Selecting the dataset:</b> Users can select a dataset in the table. The Dataset info shows the information of the dataset of interest and users can load the data by clicking the Load dataset button.

<br><br><b>Upload data:</b> Under the Upload data panel users can upload their transcriptomics and proteomics data to the platform. The platform requires 3 data files as listed below: a data file containing counts/expression (counts.csv), a sample information file (samples.csv) and a file specifying the statistical comparisons as contrasts (contrasts.csv). It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, row names and column names match for all files. On the left side of the panel, users need to provide a unique name and brief description for the dataset while uploading. N.B. Users can now create contrasts from the platform itself, so the contrasts.csv file is optional.

<br><br>
<ol>
<li>counts.csv: Count/expression file with gene on rows, samples as columns.
<li>samples.csv: Samples file with samples on rows, phenotypes as columns.
<li>contrasts.csv: Contrast file with conditions on rows, contrasts as columns.
</ol>

<br><br><br>
<center><iframe width="560" height="315" src="https://www.youtube.com/embed/elwT6ztt3Fo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe><center>

')
    
    ##-----------------------------------------------------------------------------
    ## Show current dataset on each page
    ##-----------------------------------------------------------------------------
    curDataSet <- shiny::reactive({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        ##HTML("<b>dataset :</b>",ngs$name,"")
        ##HTML("<h3>",ngs$name,"</h3>")
        dname <- gsub("^.*/|[.]pgx","",ngs$name)
        shiny::HTML("<div class='current-data'>",dname,"</div>")
    })
    output$current_dataset <- shiny::renderText({ curDataSet() })

    ##-----------------------------------------------------------------------------
    ## User interface
    ##-----------------------------------------------------------------------------
    downloadButton2 <- function (outputId, label = "Download", class = NULL, ...) {
        aTag <- shiny::tags$a(id = outputId,
                       class = paste("btn btn-default shiny-download-link", class),
                       href = "", target = "_blank", download = NA, 
                       shiny::icon("file-csv"), label, ...)
    }
    
    output$inputsUI <- shiny::renderUI({        

        delete_button <- NULL
        if(enable_delete) {
            delete_button <- shinyBS::tipify( shiny::actionButton(
                ns("deletebutton"), label=NULL, icon=icon("trash"),
                style='padding:2px 1px 1px 1px; font-size:140%; color: #B22222; width:30px;'
            ),"Delete the selected dataset.", placement="bottom")
        }

        ui <- shiny::tagList(
            shinyalert::useShinyalert(),  # Set up shinyalert
            shiny::p(shiny::strong("Dataset info:")),
            shiny::div(shiny::htmlOutput(ns("dataset_info")), id="datainfo"),
            shiny::br(),
            shiny::conditionalPanel(
                "output.rowselected != 0", ns=ns,
                shinyBS::tipify( shiny::actionButton(ns("loadbutton"),label="Load",class="load-button"),
                   "Click to load the selected dataset.", placement="bottom"),
                shinyBS::tipify( shiny::downloadButton(
                    ns("downloadpgx"), label=NULL, ## icon=icon("download"),
                    style='padding:2px 1px 1px 1px; font-size:140%; width:30px;'
                ),"Download PGX file (binary).", placement="bottom"),
                shinyBS::tipify( downloadButton2(
                    ns("downloadzip"), label=NULL, icon=icon("file-csv"),
                    style='padding:2px 1px 1px 1px; font-size:140%; width:30px;'
                ),"Download CSV files (counts.csv, samples.csv, contrasts.csv).",
                placement="bottom"),
                delete_button
            ),
            shiny::br(),shiny::br(),
            shinyBS::tipify( shiny::actionLink(ns("showfilter"), "show filters",
                                               icon=icon("cog", lib = "glyphicon")),
                   "Show dataset filters.", placement="top"),
            shiny::br(),br(),
            shiny::conditionalPanel(
                "input.showfilter % 2 == 1", ns=ns,
                shiny::uiOutput(ns("dataset_filter"))
            )
        )
        ui
    })
    shiny::outputOptions(output, "inputsUI", suspendWhenHidden=FALSE)
    output$rowselected <- shiny::reactive({
        !is.null(selectedPGX()) && length(selectedPGX())>0
    })
    shiny::outputOptions(output, "rowselected", suspendWhenHidden=FALSE)
    
    output$dataset_info <- shiny::renderText({
        sec <- currentSection()
        inf <- selectedDataSetInfo()
        inf["conditions"] <- gsub("[,]"," ",inf["conditions"])
        inf <- sapply(inf, function(s) substring(s,1,500)) 
        if(sec=="upload-data") {
            shiny::HTML(paste("<p>Please upload dataset<br>"))
        } else if(length(inf)==0) {
            shiny::HTML(paste("<p>Please select a dataset<br>"))
        } else {
            shiny::HTML(paste("<p><b>",names(inf),"</b>:", inf,""))
        }
    })

    output$dataset_filter <- shiny::renderUI({
        df <- getPGXINFO()
        datatypes <- sort(setdiff(df$datatype,c(NA,"")))
        organisms <- sort(setdiff(df$organism,c(NA,"")))        
        shiny::tagList(
            shiny::checkboxGroupInput(ns("flt_datatype"),"datatype", choices = datatypes),
            shiny::checkboxGroupInput(ns("flt_organism"),"organism", choices = organisms)
            ##checkboxGroupInput("flt_conditions","conditions",
            ##choices=c("treatment","sex","activated"))
        )

    })
    
    ##-----------------------------------------------------------------------------
    ## READ initial PGX file info
    ##-----------------------------------------------------------------------------
    ##PGXINFO <- pgx.updateInfoFile(PGX.DIR, file="datasets-info.csv", 
    ##                           force=FALSE, verbose=TRUE )
    
    ##PGXINFO  <- shiny::reactiveVal(NULL)    
    ##info <- pgx.scanInfoFile(PGX.DIR, file="datasets-info.csv", verbose=TRUE )
    ##cat("[LoadingBoard] dim.info = ",dim(info),"\n")
    ##PGXINFO(info)
    
    ## reactive value for updating table
    reload_pgxdir <- shiny::reactiveVal(0)
    
    getPGXDIR <- shiny::reactive({
        reload_pgxdir()  ## force reload
        if(!auth$logged()) return(NULL)

        email="../me@company.com"
        email <- auth$email()
        email <- gsub(".*\\/","",email)
        pdir <- pgx_dir  ## from module input
        ##USERDIR=FALSE
        if(enable_userdir) {
            pdir <- paste0(pdir,"/",email)
            if(!is.null(email) && !is.na(email) && email!="") pdir <- paste0(pdir,'/')
            if(!dir.exists(pdir)) {
                dbg("[LoadingBoard:getPGXDIR] userdir does not exists. creating pdir = ",pdir)
                dir.create(pdir)
                dbg("[LoadingBoard:getPGXDIR] copy example pgx")                
                file.copy(file.path(pgx_dir,"example-data.pgx"),pdir)
            }
        }
        pdir
    })
    
    getPGXINFO <- shiny::reactive({
        req(auth)
        if(!auth$logged()) {
            dbg("[LoadingBoard:getPGXINFO] user not logged in!")            
            return(NULL)
        }
        info <- NULL
        pdir <- getPGXDIR()
        info <- pgx.scanInfoFile(pdir, file="datasets-info.csv", verbose=TRUE )
        dbg("[LoadingBoard:getPGXINFO] pdir = ",pdir)
        dbg("[LoadingBoard:getPGXINFO] dim.info = ",dim(info))            
        if(is.null(info)) {
            aa <- rep(NA,9)
            names(aa) = c("dataset","datatype","description","nsamples",
                          "ngenes","nsets","conditions","organism","date")
            info <- data.frame(rbind(aa))[0,]        
        }        
        info
    })

    getFilteredPGXINFO <- shiny::reactive({ 
        
        ## get the filtered table of pgx datasets
        dbg("[LoadingBoard:getFilteredPGXINFO] reacted")

        req(auth)
        if(!auth$logged()) {
            dbg("[LoadingBoard:getFilteredPGXINFO] user not logged in!")            
            return(NULL)
        }
        df <- getPGXINFO()
        if(is.null(df)) return(NULL)    

        pgxdir <- getPGXDIR()
        pgxfiles = dir(pgxdir, pattern=".pgx$")
        sel <- sub("[.]pgx$","",df$dataset) %in% sub("[.]pgx$","",pgxfiles)
        df  <- df[sel,,drop=FALSE]

        ## Apply filters
        if(nrow(df)>0) {
            f1=f2=f3=rep(TRUE,nrow(df))
            notnull <- function(x) !is.null(x) && length(x)>0 && x[1]!="" && !is.na(x[1])
            if(notnull(input$flt_datatype)) f2 <- (df$datatype %in% input$flt_datatype)
            if(notnull(input$flt_organism)) f3 <- (df$organism %in% input$flt_organism)
            df <- df[which(f1 & f2 & f3),,drop=FALSE]        
            df$date <- as.Date(df$date, format='%Y-%m-%d')
            ## df$date  <- NULL
            df <- df[order(df$dataset),,drop=FALSE]   ## sort alphabetically...
            ##df <- df[order(df$date,decreasing=FALSE),]
            df <- df[order(df$date,decreasing=TRUE),]
            rownames(df) <- nrow(df):1
        }
            
        ##kk = unique(c("dataset","datatype","organism","description",colnames(df)))
        kk = unique(c("dataset","datatype","description","nsamples",
                      "ngenes","nsets","conditions","organism","date"))
        kk = intersect(kk,colnames(df))
        df = df[,kk,drop=FALSE]               
        df
    })
    
    selectedPGX <- shiny::reactive({
        ##sel <- input$pgxtable_rows_selected
        sel <- pgxtable$rows_selected()
        if(is.null(sel) || length(sel)==0) return(NULL)
        df <- getFilteredPGXINFO()
        if(is.null(df) || nrow(df)==0) return(NULL)        
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$","",pgxfile),".pgx") ## add/replace .pgx
        pgxfile
    })

    selectedDataSetInfo <- shiny::reactive({
        ##sel <- input$pgxtable_rows_selected
        sel <- pgxtable$rows_selected()
        if(is.null(sel) || length(sel)==0) return(NULL)
        df <- getFilteredPGXINFO()
        if(is.null(df) || nrow(df)==0) return(NULL)        
        unlist(lapply(df[sel,],as.character))
    })

    ##=============================================================================
    ##========================== OBSERVE/REACT ====================================
    ##=============================================================================
    pgxfile="geiger2016-arginine"

    loadPGX <- function(pgxfile) {        
        pgxfile <- paste0(sub("[.]pgx$","",pgxfile),".pgx") ## add/replace .pgx         
        pgxdir <- getPGXDIR()
        pgx.path <- pgxdir[file.exists(file.path(pgxdir,pgxfile))][1]
        pgxfile1 = file.path(pgx.path,pgxfile)
        pgxfile1
        ngs <- NULL
        pgx <- NULL
        if(file.exists(pgxfile1)) {
            message("[LoadingBoard::loadPGX] loading ",pgxfile)
            shiny::withProgress(message="Loading data...", value=0.33, {            
                load(pgxfile1,verbose=0)
            })
        } else {
            cat("[LoadingBoard::loadPGX] error file not found : ",pgxfile)
            return(NULL)
        }
        if(!is.null(ngs)) return(ngs)
        if(!is.null(pgx)) return(pgx)
    }
    
    output$downloadpgx <- shiny::downloadHandler(
        ##filename = "userdata.pgx",
        filename = function() {
            selectedPGX()
        },
        content = function(file) {
            pgxfile <- selectedPGX()
            cat("[LoadingBoard::loadPGX] pgxfile = ",pgxfile,"\n")
            if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)
            ngs <- loadPGX(pgxfile)
            temp <- tempfile()
            cat("[LoadingBoard::loadPGX] temp = ",temp)            
            save(ngs, file=temp)
            file.copy(temp,file)
        }
    )

    output$downloadzip <- shiny::downloadHandler(
        ##filename = "userdata.zip",
        filename = function() {
            sub("pgx$","zip",selectedPGX())
        },
        content = function(file) {
            ## ngs <- currentPGX()  ## current dataset
            pgxfile <- selectedPGX()
            cat("[LoadingBoard::downloadZIP] pgxfile = ",pgxfile,"\n")            
            if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)
            pgxname <- sub("[.]pgx$","",pgxfile)
            ngs <- loadPGX(pgxfile)
            dir.create(tmp <- tempfile())
            tmp2 <- file.path(tmp,pgxname)
            dir.create(tmp2)

            exp.matrix <- sign(ngs$model.parameters$exp.matrix)
            exp.matrix <- contrastAsLabels(exp.matrix) ## new recommended style
            exp.matrix[is.na(exp.matrix)] <- ""
            
            write.csv(ngs$counts,  file=file.path(tmp2, "counts.csv"))
            write.csv(ngs$samples, file=file.path(tmp2, "samples.csv"))
            write.csv(exp.matrix, file=file.path(tmp2, "contrasts.csv"))
            zipfile <- tempfile(fileext = ".zip")
            zip::zip(zipfile,
                     files=paste0(pgxname,"/",c("counts.csv","samples.csv","contrasts.csv")),
                     root=tmp)
            ## zip::zip_list(zipfile)
            cat("[LoadingBoard::downloadZIP] zipfile = ",zipfile)                        
            file.copy(zipfile,file)
        }
    )
    
    shiny::observeEvent( input$deletebutton, {
                
        ##pgxfile <- currentPGX()$name
        pgxfile <- selectedPGX()
        if(is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) return(NULL)
        ##pgx.path <- PGX.DIR[file.exists(file.path(PGX.DIR,pgxfile))][1]
        pgx.path <- getPGXDIR()
        pgxfile1 = file.path(pgx.path,pgxfile)
        pgxfile1
        sel <- NULL

        deletePGX <- function() {
            if(input$confirmdelete) {
                cat(">>> deleting",pgxfile,"\n")
                pgxfile2 <- paste0(pgxfile1,"_")  ## mark as deleted
                file.rename(pgxfile1, pgxfile2)
                ##touchtable(touchtable()+1)
                reload_pgxdir(reload_pgxdir()+1)                
                if(0) {
                    this.pgx <- sub("[.]pgx$","",pgxfile)
                    all.pgx  <- sub("[.]pgx$","",PGXINFO()$dataset)
                    table.pgx <- sub("[.]pgx$","",getFilteredPGXINFO()$dataset)
                    sel <- which(table.pgx == this.pgx)
                    newpgx <- PGXINFO()[all.pgx != this.pgx,]
                    PGXINFO(newpgx)
                    DT::selectRows(proxy = DT::dataTableProxy(ns("pgxtable")), selected=sel)
                }
                    
            } else {
                cat(">>> deletion cancelled\n")
            }
        }

        not.anonymous <- !is.na(auth$name()) && auth$name()!=""         
        allow.delete <- !not.anonymous        
        message("[LoadingBoard::@deletebutton] current user = ",auth$name()," \n")
        message("[LoadingBoard::@deletebutton] allow.delete = ",allow.delete," \n")
        
        allow.delete = TRUE
        if(!allow.delete) {
            message("[LoadingBoard::@deletebutton] WARNING:: ",pgxfile,
                    " not owned by ",auth$name()," \n")
            shinyalert::shinyalert(
                            title = "Error!",
                            text = "You do not have permission to delete this dataset",
                            type = "error"
                        )
        } else {
            shinyalert::shinyalert(
                "Delete this dataset?",
                paste("Are you sure you want\nto delete '",pgxfile,"'?"),
                confirmButtonText = "Delete",
                showCancelButton = TRUE,
                callbackR = deletePGX,
                inputId = "confirmdelete")
        }

        
    })
    
    ##=================================================================================
    ##=============================== MODAL DIALOGS ===================================
    ##=================================================================================
    
    startup_count=0

    particlesjs.conf <- rjson::fromJSON(file="resources/particlesjs-config.json")

    showStartupModal <- function(once=FALSE) {    
        if(length(input$loadbutton)==0) {
            dbg("[showStartupModal] UI not ready. skipping")
            ##delay(8000, DT::selectRows(proxy = DT::dataTableProxy("pgxtable"), selected=1))
            ##delay(8000, shinyjs::click("loadbutton"))
            return(NULL)  ## UI not ready???rt
        }
        if(once && startup_count>0) return(NULL)
        
        dbg("showStartupModal: showing!\n")    
        AuthenticationUI(ns("auth"))
        
        startup_count <<- startup_count + 1    
        dbg("showStartupModal done!\n")
    }

    ## shiny::observeEvent( input$action_beer, {
    ##     dbg("buy beer button action\n")
    ##     startup_count <<- startup_count + 1    
    ##     USER$logged <- TRUE
    ##     USER$name   <- "beer buddy"
    ##     shiny::removeModal()
    ##     ##alert("Wow. Thanks buddy!")
    ##     shinyWidgets::sendSweetAlert(
    ##         session=session, title="Wow. Thanks buddy!",
    ##         text = "Free entrance for you!", type = "info")
    ##     ##Sys.sleep(4);removeModal()
    ## })

    currentSection <- shiny::reactive({
        cdata <- session$clientData
        sub("section-","",cdata[["url_hash"]])
    })

    ##=================================================================================
    ##========================= USER AUTHENTICATION ===================================
    ##=================================================================================

    ## USER <- shiny::reactiveValues( logged = FALSE, name="anonymous")
    
    ##================================================================================
    ##====================== INPUT DATA REACTIVE OBJECT ==============================
    ##================================================================================

    currentPGX <- shiny::reactiveVal(NULL)
    
    ##inputData <- shiny::eventReactive( reload(), {
    inputData <- shiny::reactive({
        ## This is wrapper function that loads the ngs object. It also
        ## checks if the user is logged in.
        ##  
        dbg("[LoadingBoard::inputData] !!! reacted !!!\n")
        dbg("[LoadingBoard::inputData] authentication=",authentication,"\n")
        dbg("[LoadingBoard::inputData] auth$logged=",auth$logged(),"\n")
        
        if(!auth$logged()) {
            currentPGX(NULL)  ## set NULL
            return(NULL)
        }
        
        pgx <- currentPGX()            
        dbg("[LoadingBoard::inputData] pgx$name = ",pgx$name,"\n")        
        return(pgx)
    })

    shiny::observeEvent( input$loadbutton, {

        ## Observe button press
        btn <- shiny::isolate(input$loadbutton)
        pgxfile = NULL
        pgxfile = shiny::isolate(selectedPGX())

        if(is.na(pgxfile) || is.null(pgxfile) || pgxfile=="" || length(pgxfile)==0) {
            return(NULL)
        }

        if(!is.null(btn) && btn!=0) {
            ## show loading pop-up
            pgx.showCartoonModal()
        }
        
        message("[LoadingBoard::<loadbutton>] loading pgxfile = ",pgxfile)
        pgx <- loadPGX(pgxfile)

        if(is.null(pgx)) {
            cat("[LoadingBoard::<loadbutton>] ERROR file not found : ",pgxfile,"\n")
            beepr::beep(10)
            if(loadedDataset()) shiny::removeModal()
            return(NULL)
        }
        
        ##----------------- update input
        message("[LoadingBoard::<loadbutton>] initializing PGX object")
        pgx <- pgx.initialize(pgx)
        message("[LoadingBoard::<loadbutton>] initialization done!")        
        if(is.null(pgx)) {
            cat("[LoadingBoard::<loadbutton>] ERROR in object initialization\n")
            beepr::beep(10)
            shiny::showNotification("ERROR in object initialization!\n")
            if(loadedDataset()) shiny::removeModal()
            return(NULL)
        }
        if(is.null(pgx$name)) pgx$name <- sub("[.]pgx$","",pgxfile)

        ##----------------- remove modal??
        if(startup_count>0) {
            Sys.sleep(4)
            if(loadedDataset()) shiny::removeModal()
        }
        
        loadedDataset(TRUE)
        currentPGX(pgx)
    })
    ##}, ignoreNULL=FALSE )
    ##}, ignoreNULL=TRUE )


    ## ================================================================================
    ## =============================== VALUE BOXES UI =================================
    ## ================================================================================


    ## shinyWidgets::useShinydashboard()
    vbox <- function(value, label) {
        shinydashboard::box(
            shiny::h1(value, style="font-weight: 800; color: white; padding: 16px 0 0 0; margin: 0 0 0 20px;"),
            shiny::h5(label, style="margin: 0 0 0 20px; font-weight: 400; color: white; padding-bottom: 25px"),
            ##h1(value, style="font-weight: 800; padding: 16px 0 0 0; margin: 0 0 0 20px;"),
            ##h5(label, style="margin: 0 0 0 20px; font-weight: 400; padding-bottom: 25px"),
            width="100%", class="vbox")
    }

    output$valuebox1 <- shiny::renderUI({
        pgx <- getFilteredPGXINFO()
        shiny::req(pgx)
        ndatasets = "..."
        ndatasets <- nrow(pgx)
        vbox( ndatasets, "data sets")     
    })

    output$valuebox2 <- shiny::renderUI({
        pgx <- getFilteredPGXINFO()
        shiny::req(pgx)
        ##dbg("valuebox2:: pgx$nsamples=",pgx$nsamples)
        nsamples <- sum(as.integer(pgx$nsamples),na.rm=TRUE)
        vbox( nsamples, "number of samples")     
    })

    output$valuebox3 <- shiny::renderUI({
        pgx <- getFilteredPGXINFO()
        shiny::req(pgx)
        ##dbg("valuebox3:: pgx$nsamples=",pgx$nsamples)
        nvalues <- sum(as.integer(pgx$nsamples) * (as.integer(pgx$ngenes)
            + as.integer(pgx$nsets)),na.rm=TRUE)
        nvalues1 <- format(nvalues, nsmall=, big.mark=" ")  # 1,000.6
        vbox( nvalues1, "data points") 
    })

    output$valueboxes_UI <- shiny::renderUI({
        shiny::fillRow(
            height=115,
            shiny::uiOutput(ns("valuebox1")),
            shiny::uiOutput(ns("valuebox2")), 
            shiny::uiOutput(ns("valuebox3"))
        )
    })

    ##================================================================================
    ## Data sets
    ##================================================================================
    
    ## reactive value for updating table
    touchtable <- shiny::reactiveVal(0)

    ##split=" ";n=5
    andothers <- function(s, split=" ", n=8) {
        if(is.na(s)) return("")
        s <- sub("^[ ]*","",s)
        s <- sub("[ ]+"," ",s)
        s1 <- strsplit(s, split=split)[[1]]
        if(length(s1)<=n) return(s)
        n2 <- setdiff(length(s1),n)
        paste(paste(head(s1,n), collapse=" "),"(+",n2,"others)")
    }

    pgxTable.RENDER <- shiny::reactive({

        if(SHOWSPLASH) showStartupModal(once=TRUE)
        
        dbg("[pgxTable.RENDER] reacted")

        ##touchtable()  ## explicit reactive on this
        reload_pgxdir()
        
        df <- getFilteredPGXINFO()
        shiny::req(df)
        dbg("[pgxTable.RENDER] dim(df)=",dim(df))

        updateTab <- function() {
            ## NEED RETHINK!!! DOES NOT WORK!!!
            dbg("[pgxTable.RENDER] updating tab to Upload Data")            
            updateTabsetPanel(session, ns("tabs"), selected = "Upload data")
            updateTabsetPanel(session, "load-tabs", selected = "Upload data")
        }

        if(FALSE && nrow(df)==0 && auth$logged()) {
            ## NEED RETHINK. Sometimes pops up at login...
            dbg("[pgxTable.RENDER] nrow(df) = ",nrow(df))
            dbg("[pgxTable.RENDER] auth$logged() = ",auth$logged())            
            shinyalert::shinyalert(
                            title = "Your playground looks empty...",
                            text = "Please start by uploading some data!",
                            type = "warning",
                            callbackR = updateTab
                        )
        }
        
        df$dataset  <- gsub("[.]pgx$"," ",df$dataset)
        df$conditions  <- gsub("[,]"," ",df$conditions)
        df$conditions  <- sapply(as.character(df$conditions), andothers, split=" ", n=5)
        df$description <- shortstring(as.character(df$description),200)
        df$nsets <- NULL

        target1 <- grep("date",colnames(df))
        target2 <- grep("description",colnames(df))
        
        DT::datatable(df,
                      class = 'compact cell-border stripe hover',
                      rownames=TRUE,
                      extensions = c('Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      ## filter = "top",
                      fillContainer = TRUE,
                      options=list(
                          ##dom = 'Blfrtip',
                          dom = 'frti',
                          ##columnDefs = list(list(searchable = FALSE, targets = 1)),
                          pageLength = 1000, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = FALSE,
                          ##scrollY =400, ## scroller=TRUE,
                          scrollY = '100vh', ## scroller=TRUE,
                          deferRender=TRUE,
                          autoWidth = TRUE,
                          columnDefs = list(
                              list(width='50px', targets=target1),
                              list(width='260px', targets=target2)
                          )
                      )  ## end of options.list 
                      )  %>%
            DT::formatStyle(0, target='row', fontSize='11.5px', lineHeight='95%')

    })

    pgxtable_text = "This table contains a general information about all available datasets within the platform. For each dataset, it reports a brief description as well as the total number of samples, genes, gene sets (or pathways), corresponding phenotypes and the creation date."

    pgxtable <- shiny::callModule(
        tableModule, id = "pgxtable",
        func = pgxTable.RENDER,
        title = "Datasets",
        height = 640, width = c('100%',1600),
    )

    output$pgxtable_UI <- shiny::renderUI({    
        shiny::fillCol(
            height = 750,
            flex = c(NA,1),
            shiny::uiOutput(ns("valueboxes_UI")),
            shiny::fillRow(
                flex = c(1,0.1,4.5),
                shiny::wellPanel(
                    shiny::uiOutput(ns("inputsUI"))
                ),
                shiny::br(), 
                tableWidget(ns("pgxtable"))
            )
        )        
    })
    shiny::outputOptions(output, "pgxtable_UI", suspendWhenHidden=FALSE) ## important!

    ##================================================================================
    ## Upload new data
    ##================================================================================
    if(enable_upload) {
        dbg("[LoadingBoard] upload enabled")
        
        output$upload_UI <- shiny::renderUI({
            dbg("[LoadingBoard] output$upload_UI::renderUI")
            UploadModuleUI(ns("upload_panel"))
        })

        dbg("[LoadingBoard] initializing UploadModule")        
        uploaded_pgx <- UploadModuleServer(
            id = "upload_panel",
            FILES = FILES,
            pgx.dirRT = shiny::reactive(getPGXDIR()),
            height = 720,
            ## limits = c(samples=20, comparisons=20, genes=8000),
            limits = limits
        )
        
        shiny::observeEvent( uploaded_pgx(), {
            
            dbg("[observe::uploaded_pgx] uploaded PGX detected!")
            pgx <- uploaded_pgx()
            
            dbg("[observe::uploaded_pgx] initializing PGX object")
            pgx <- pgx.initialize(pgx)
            
            ## update CurrentPGX
            currentPGX(pgx)
            DT::selectRows(proxy = DT::dataTableProxy(ns("pgxtable")), selected=NULL)
            
            savedata_button <- NULL
            if(enable_save) {                
                dbg("[LoadingBoard] observeEvent:savedata reacted")        
                ## -------------- save PGX file/object ---------------
                ##pgx <- currentPGX()
                pgxname <- sub("[.]pgx$","",pgx$name)
                pgxname <- paste0(gsub("[ \\/]","_",pgxname),".pgx")
                pgxdir  <- getPGXDIR()
                fn  <- file.path(pgxdir,pgxname)
                
                ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ## Note: Currently we use 'ngs' as object name but want to go
                ## towards 'pgx' as standard name. Actually saving as RDS
                ## should be better.
                ngs=pgx
                save(ngs, file=fn)
                remove(ngs)
                
                if(1) {
                    message("[LoadingBoard::@savedata] updating PGXINFO")                    
                    pgx.initDatasetFolder(pgxdir, force=FALSE, verbose=TRUE)
                } else {
                    message("[LoadingBoard::@savedata] updating PGXINFO file")
                    old.info <- getPGXINFO()
                    new.info <- pgx.updateInfoPGX(old.info, pgx, remove.old=TRUE)
                    pgxdir <- getPGXDIR()[1] ## first folder!!
                    info.file <- file.path(pgxdir, "datasets-info.csv")  
                    Sys.chmod(info.file, mode="0666")
                    write.csv(new.info, file=info.file)
                    message("[LoadingBoard::@savedata] saved PGXINFO file!")                
                }
                
                ##PGXINFO(new.info)
                ##touchtable(touchtable()+1)  ## refresh table
                reload_pgxdir(reload_pgxdir()+1)
            }
            
            ## shiny::removeModal()
            msg1 <- "<b>Ready!</b>"
            ##beepr::beep(sample(c(3,4,5,6,8),1))  ## music!!
            beepr::beep(10)  ## short beep
            
            if(enable_save) {
                msg1 <- "<b>Ready!</b><br>Your data is ready and has been saved in your library. You can now start exploring your data."
            } else {
                msg1 <- "<b>Ready!</b><br>Your data is ready. You can now start exploring your data."
            }
            loadedDataset(TRUE)
            shiny::showModal(
                       shiny::modalDialog(
                                  shiny::HTML(msg1),
                                  title = NULL,
                                  size = "s",
                                  footer = shiny::tagList(
                                                      ## savedata_button,
                                                      ## shiny::actionButton(ns("sharedata"), "Share with others", icon=icon("share-alt")),
                                                      shiny::modalButton("Start!")
                                                  )
                              ))
            shiny::updateTabsetPanel(session, "tabs",  selected = "Datasets")
        })

    }

    ##---------------------------------------------------------------
    ##----------------- modules for Forum ---------------------------
    ##---------------------------------------------------------------
    
    output$forum <- shiny::renderUI({
        parenturl <- paste0(session$clientData$url_protocol,
                            "//",session$clientData$url_hostname,
                            ":",session$clientData$url_port,
                            session$clientData$url_pathname)
        ## parenturl <- gsub("localhost","127.0.0.1",parenturl)
        parenturl <- URLencode(parenturl, TRUE)
        cat("[LoadingBoard:forum] parenturl =",parenturl,"\n")
        src = paste0('https://groups.google.com/forum/embed/?place=forum/omicsplayground',
                     '&showsearch=true&showpopout=true&parenturl=',parenturl)
        cat("src = ",src,"\n")
        shiny::tags$iframe(id="forum_embed", src=src, height=600, width='100%',
                    ##seamless="seamless",
                    frameborder='no')
        ##HTML(src)
    })
         
    output$tweet <- shiny::renderUI({
        ## NOT WORKING YET...
        shiny::tags$a(class="twitter-timeline",
               href="https://twitter.com/bigomics?ref_src=twsrc%5Etfw")
        ##shiny::tags$script('twttr.widgets.load(document.getElementById("tweet"));')
    })
            
    output$forum_UI <- shiny::renderUI({
        shiny::fillCol(
            height = 550,
            shiny::fillRow(
                flex=c(4,0),
                shiny::htmlOutput(ns("forum"))
                ##uiOutput("tweet")
            )
        )
    })
    
    ##------------------------------------------------
    ## Board return object
    ##------------------------------------------------
    res <- list(
        loaded = loadedDataset,
        inputData = inputData,
        auth = auth
        ##inputData = currentPGX,
        ##usermode = shiny::reactive({ USERMODE() })
    )
    return(res)
}
