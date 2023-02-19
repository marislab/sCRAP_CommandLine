
suppressPackageStartupMessages(library(optparse))

    option_list <- list(
      make_option(c("-p", "--peptide"),default = "", type = "character",
                  help = "(Required) Peptide query - typically about 9 peptides {character string}"),
      make_option(c("-l", "--hla"),default = "",  type = "character",
                  help = "(Required) hla allele in format: HLA-{allele},  eg.'HLA-A0201' 'HLA-A0301' "),
      make_option(c("-s", "--hotspots"), type = "character",default = "4,5,6",
                  help = "hotspots locations, numeric value 1-9 seperated by commas - {default: 4,5,6}"),
      make_option(c("-o", "--outdir"), type = "character",
                  help = "Output directory path", default = "."),
      make_option(c("-f", "--filename"), default = "", type = "character",
                  help = "Output file basename"),
      make_option(c("-d", "--datadir"), default = "data/", type = "character",
                  help = "data directory path"),
      make_option(c("-m", "--hladir"), default = "HLA/", type = "character",
                 help = "hla data directory path"),
      make_option(c("-t", "--filetype"), default = "csv", type = "character",
                  help = "select output file format (csv/tsv) {default: csv}")
    )
    
    CHECKS = TRUE
    
    
    # parse the parameters
    CHECKS=tryCatch(expr={
      parse_args(OptionParser(option_list = option_list,add_help_option = T),print_help_and_exit = T); 
      TRUE
    }, error=function(e){
      cat("\nThere was an error with your input arguments")
      return(F)
      }
    )# tryCatch({


    if(CHECKS)
    {
      opt <- parse_args(OptionParser(option_list = option_list,add_help_option = T),print_help_and_exit = T)
  
      
      hla_options = c("HLA-A0101","HLA-A0201","HLA-A0202","HLA-A0203","HLA-A0205","HLA-A0206","HLA-A0207","HLA-A0211","HLA-A0212","HLA-A0216","HLA-A0217","HLA-A0219","HLA-A0250","HLA-A0301","HLA-A0302","HLA-A0319","HLA-A1101","HLA-A2301","HLA-A2402","HLA-A2403","HLA-A2501","HLA-A2601","HLA-A2602","HLA-A2603","HLA-A2902","HLA-A3001","HLA-A3002","HLA-A3101","HLA-A3201","HLA-A3207","HLA-A3215","HLA-A3301","HLA-A6601","HLA-A6801","HLA-A6802","HLA-A6823","HLA-A6901","HLA-A8001","HLA-B0702","HLA-B0801","HLA-B0802","HLA-B0803","HLA-B1401","HLA-B1402","HLA-B1501","HLA-B1502","HLA-B1503","HLA-B1509","HLA-B1517","HLA-B1801","HLA-B2705","HLA-B2720","HLA-B3501","HLA-B3503","HLA-B3701","HLA-B3801","HLA-B3901","HLA-B4001","HLA-B4002","HLA-B4013","HLA-B4201","HLA-B4402","HLA-B4403","HLA-B4501","HLA-B4506","HLA-B4601","HLA-B4801","HLA-B5101","HLA-B5301","HLA-B5401","HLA-B5701","HLA-B5703","HLA-B5801","HLA-B5802","HLA-B7301","HLA-B8101","HLA-B8301","HLA-C0303","HLA-C0401","HLA-C0501","HLA-C0602","HLA-C0701","HLA-C0702","HLA-C0802","HLA-C1203","HLA-C1402","HLA-C1502")
      
      
      input_PepIn <- opt$peptide
      input_HLA <- opt$hla
      input_hotspots <- unlist(strsplit(opt$hotspots,split=","))
      input_outdir <- opt$outdir
      input_filename <- opt$filename
      input_datadir <- opt$datadir
      input_hladir <- opt$hladir
      input_filetype <- tolower(opt$filetype)
    
      if(nchar(input_PepIn) <= 3)
      {
        cat("\nError: Enter a peptide sequence\n")
        CHECKS = F
      }
  
      if(length(which(input_HLA %in% hla_options)) < 1)
      {
        cat("\nError: Enter a valid HLA sequence with the '-l' option (See below).\nValid options inlcude:\n")
        cat(paste("  ",paste(hla_options[1:10],collapse = " "),"\n"))
        cat(paste("  ",paste(hla_options[11:20],collapse = " "),"\n"))
        cat(paste("  ",paste(hla_options[21:30],collapse = " "),"\n"))
        cat(paste("  ",paste(hla_options[31:40],collapse = " "),"\n"))
        cat(paste("  ",paste(hla_options[41:50],collapse = " "),"\n"))
        cat(paste("  ",paste(hla_options[51:60],collapse = " "),"\n"))
        cat(paste("  ",paste(hla_options[61:70],collapse = " "),"\n"))
        cat(paste("  ",paste(hla_options[71:80],collapse = " "),"\n"))
        cat(paste("  ",paste(hla_options[80:87],collapse = " "),"\n"))
        
        CHECKS = F
      }
  
      if(input_hotspots[1] == "0")
      {
        input_hotspots = c()
      }
      
      if(length(input_hotspots) >= 1)
      {
        input_hotspots = as.numeric(input_hotspots)[which(as.numeric(input_hotspots) > 0 & as.numeric(input_hotspots) <= 9) ]
      }
      
      if(length(input_hotspots) >= 1)
      {
        input_hotspots = as.numeric(input_hotspots)[which(as.numeric(input_hotspots) > 0 & as.numeric(input_hotspots) <= 9) ]
      }
      
      if(input_filetype != "csv" & input_filetype != "tsv")
      {
        input_filetype = "csv"
      }
      
      DATA_File_Check = c("allelelist.RDS","gtex_annot.RDS","genekey.RDS","gtex_p0.RDS","gtex_p1.RDS","gtex_p2.RDS","gtex_p3.RDS","gtex_p4.RDS","gtex_p5.RDS","gtex_p6.RDS","gtex_p7.RDS","gtex_p8.RDS","gtex_p9.RDS","combined_normal_ligands.RDS")
      
      HLAfiles <- list.files(path=input_hladir, pattern = "*.RDS")
      Datafiles <- list.files(path=input_datadir, pattern = "*.RDS")
      
      if(!all(DATA_File_Check %in% Datafiles))
      {
        CHECKS = F
        
        cat("\nError: Please enter a valid data directory and ensure all the required files are present:\n")
        cat(paste("  ",paste(DATA_File_Check[1:3],collapse = " "),"\n"))
        cat(paste("  ",paste(DATA_File_Check[4:6],collapse = " "),"\n"))
        cat(paste("  ",paste(DATA_File_Check[7:9],collapse = " "),"\n"))
        cat(paste("  ",paste(DATA_File_Check[10:12],collapse = " "),"\n"))
        cat(paste("  ",paste(DATA_File_Check[13:14],collapse = " "),"\n"))
        
      }
  
      if(nchar(input_filename) == 0  )
      {
        out_file = paste("Result_",input_PepIn,"_",input_HLA,".",input_filetype,sep="")
      }
      if(nchar(input_filename) >= 1  )
      {
        out_file = paste(input_filename,".",input_filetype,sep="")
      }
      
      Output_Path = paste(input_outdir,"/",out_file,sep="")
      
    } #if(CHECKS)  
    
    if(!CHECKS)
    {
      cat(paste("\nUsage: Rscript sCRAP.R -p <peptide> -l <HLA-allele> {arguments}\n\nuse:\n\n  Rscript sCRAP.R -h \n\nfor more help.\n\n"))
    } #if(!CHECKS)
    
      
if(CHECKS){
  
  
  HLAfiles <- list.files(path=input_hladir, pattern = "*.RDS")
  
  cat(" Loading Data...\n")
  
  gtex_annot <- readRDS(paste(input_datadir,"/gtex_annot.RDS",sep=""))
  genekey <- readRDS(paste(input_datadir,"/genekey.RDS",sep=""))
  gtex_p0 <- readRDS(paste(input_datadir,"/gtex_p0.RDS",sep=""))
  gtex_p1 <- readRDS(paste(input_datadir,"/gtex_p1.RDS",sep=""))
  gtex_p2 <- readRDS(paste(input_datadir,"/gtex_p2.RDS",sep=""))
  gtex_p3 <- readRDS(paste(input_datadir,"/gtex_p3.RDS",sep=""))
  gtex_p4 <- readRDS(paste(input_datadir,"/gtex_p4.RDS",sep=""))
  gtex_p5 <- readRDS(paste(input_datadir,"/gtex_p5.RDS",sep=""))
  gtex_p6 <- readRDS(paste(input_datadir,"/gtex_p6.RDS",sep=""))
  gtex_p7 <- readRDS(paste(input_datadir,"/gtex_p7.RDS",sep=""))
  gtex_p8 <- readRDS(paste(input_datadir,"/gtex_p8.RDS",sep=""))
  gtex_p9 <- readRDS(paste(input_datadir,"/gtex_p9.RDS",sep=""))
  
  gtex <- rbind(gtex_p0,gtex_p1,gtex_p2,gtex_p3,gtex_p4,gtex_p5,gtex_p6,gtex_p7,gtex_p8,gtex_p9)
  
  combined <- readRDS(paste(input_datadir,"/combined_normal_ligands.RDS",sep=""))
  

    G=list()
    G[[1]] = c("A","G") #Short Chains
    G[[2]] = c("K","R","H") #Basic
    G[[3]] = c("N","Q") #Polar Uncharged
    G[[4]] = c("D","E") #Acidic
    G[[5]] = c("C", "M") #Contains S
    G[[6]] = c("F","Y","W","H") #Aromatic
    G[[7]] = c("S","T","Y") #Alcohol-Hydroxyl Group
    G[[8]] = c("I","L","V","M","A") #Aliphatic
    G[[9]] = c("P") #Weird - I mean proline
    
    Polarity = list()
    Polarity[[1]] = c("D","E","N","Q","R","K","H","Y","C","S","T") #Polar
    Polarity[[2]] = c("G","A","F","W","P","I","L","V","M") #Non Polar
    
    Charge = list()
    Charge[[1]] = c("D","E") #Negative
    Charge[[2]] = c("K","R","H") #Positive
    
    
    HLA_RDS_Files <- system(paste("ls HLA/", input_HLA, ".p*.RDS", sep=""), intern=T)
    
    combinedpeps <- data.frame()
    for(I in 1:length(HLA_RDS_Files))
    {
      combinedpeps <- rbind(combinedpeps,data.frame(readRDS(HLA_RDS_Files[I])))
    }
    
    
    peptable <- combinedpeps[,c(3,11,13,15)]
    peptable$Identity <- as.character(gsub("_HUMAN","", as.character(peptable$Identity)))
    
    
    #Function to Compare Peptides
    Get_Peptide_Scores = function(QUERY,PEPTIDE)
    {
      QUERY.array = unlist(strsplit(as.character(QUERY),split=""))
      PEPTIDE.array = unlist(strsplit(as.character(PEPTIDE),split=""))
      
      Peptide_Score = -100
      if(
        length(QUERY.array) == length(PEPTIDE.array) #&
      )#if
      {#Start function
        TEST.Pair=rbind(QUERY.array,PEPTIDE.array)
        
        #Equivalent
        EQU_score = apply(TEST.Pair,MARGIN = 2,function(x) {Score=0;Criteria = x[1]==x[2]; if(Criteria){Score=3};return(Score)})
        #GROUP
        GROUP_score = apply(TEST.Pair,MARGIN=2,function(X) {Score=0;Criteria=any(unlist(lapply(G,function(G) all(X %in% G))));if(Criteria){Score=2};return(Score)} )
        #Polarity
        Polarity_score = apply(TEST.Pair,MARGIN=2,function(X) {Score=0;Criteria=any(unlist(lapply(Polarity,function(G) all(X %in% G))));if(!Criteria){Score=-2};return(Score)})
        #Charge
        Charge_score = apply(TEST.Pair,MARGIN=2,function(X) {Score=0;Criteria=all(unlist(lapply(Charge, function(x) any(x %in% G))));if(!Criteria){Score=-1};return(Score)})
        
        #Double Scores for Selected Hotspots.
        
        HOTSPOTS <- as.numeric(unlist(input_hotspots))
        
        EQU_score[HOTSPOTS] = EQU_score[HOTSPOTS]*2 
        GROUP_score[HOTSPOTS] = GROUP_score[HOTSPOTS]*2 
        Polarity_score[HOTSPOTS] = Polarity_score[HOTSPOTS]*2 
        Charge_score[HOTSPOTS] = Charge_score[HOTSPOTS]*2 
        
        #Sets scores at Positions 2 and last position to 0.
        EQU_score[c(2,length(EQU_score))] = 0
        GROUP_score[c(2,length(GROUP_score))] = 0
        Charge_score[c(2,length(Charge_score))] = 0
        
        #Final Match Score
        Peptide_Score = sum(c(EQU_score,GROUP_score,Polarity_score))
        
      }#{Start function
      return(Peptide_Score)
    }
    
    #normal tissue RPKM
    normal <- function(peptide){
      for (i in 1:nrow(peptide)){
        if (!is.na(g <- match(peptide$Gene[i], genekey$name))){
          g <- match(peptide$Gene[i], genekey$name)
          peptide$Gene[i] <- genekey$gene[g]
        }
        if (!is.na(match(peptide$Gene[i], rownames(gtex)))){
          r <- match(peptide$Gene[i], rownames(gtex))
          peptide$max_norm[i] <- max(gtex[r,])
          
          tiss <- match(colnames(gtex)[which.max(gtex[r,])[[1]]], gtex_annot$SAMPID)
          peptide$max_tissue[i] <- as.character(gtex_annot$SMTS[as.numeric(tiss)])
        }
        else {
          peptide$max_norm[i] <- NA
          peptide$max_tissue[i] <- NA
        }
        

      } #for (i in 1:nrow(peptide)) (end)
      return(peptide)
    } #normal <- function(peptide) (end)
    
    # compare to  normal ligandome
    ligandome <- function(peptide){
      
      
      peptide$Normal_Ligandome <- "F"
      peptide$MHC_class <- ""
      peptide$HLA <- ""
      peptide$Organism <- ""
      for (i in 1:nrow(peptide)){
        if (!is.na(match(peptide$Peptide[i], combined$search_hit))){
          r <- match(peptide$Peptide[i], combined$search_hit)
          peptide$Normal_Ligandome[i] <- "T"
          peptide$MHC_class[i] <- as.character(combined$MHCClass[r])
          peptide$HLA[i] <- as.character(combined$top_allele[r])
          peptide$Organism[i] <- as.character(combined$Organism[r])
        }
        else{
          peptide$Normal_Ligandome[i] <- "F"
        }

      } #for (i in 1:nrow(peptide)) (end)
      
      return(peptide)
    } #ligandome <- function(peptide) (end)
    
    
    
    {
      
      PEPTIDE_INPUT = input_PepIn
      cat(" Filtering table: this can take several minutes\n")
      
      peptable.peptide.subset = peptable$peptide[which(sapply(as.character(peptable$peptide),nchar) == nchar(PEPTIDE_INPUT) &  peptable$BindLevel=="<=SB")]
      
      peptable.peptide.subset = peptable.peptide.subset[which((substr(peptable.peptide.subset, 1,1) == substr(PEPTIDE_INPUT, 1,1)) | (substr(peptable.peptide.subset, 3,3) == substr(PEPTIDE_INPUT, 3,3)) | (substr(peptable.peptide.subset, 4,4) == substr(PEPTIDE_INPUT, 4,4)) | (substr(peptable.peptide.subset, 5,5) == substr(PEPTIDE_INPUT, 5,5)) | (substr(peptable.peptide.subset, 6,6) == substr(PEPTIDE_INPUT, 6,6)) | (substr(peptable.peptide.subset, 7,7) == substr(PEPTIDE_INPUT, 7,7)) | (substr(peptable.peptide.subset, 8,8) == substr(PEPTIDE_INPUT, 8,8)))]
      
      cat(paste(" Calculating", length(peptable.peptide.subset)," scores. This may take a few minutes\n"))
      
      quarter_count = floor(length(peptable.peptide.subset)/4)
      mod_count = length(peptable.peptide.subset) %% 4
      
      TEST_Scores.filter_1 = sapply(peptable.peptide.subset[1:quarter_count],function(X) as.numeric(Get_Peptide_Scores(PEPTIDE_INPUT, X)))
      
      cat(paste(" Calculating", (length(peptable.peptide.subset)-quarter_count)," scores. This may take a few minutes\n"))
      
      TEST_Scores.filter_2 = sapply(peptable.peptide.subset[(quarter_count+1):(quarter_count*2)],function(X) as.numeric(Get_Peptide_Scores(PEPTIDE_INPUT, X)))
      
      cat(paste(" Calculating", (length(peptable.peptide.subset)-quarter_count*2)," scores. This may take a few minutes\n"))
      TEST_Scores.filter_3 = sapply(peptable.peptide.subset[(quarter_count*2 + 1):(quarter_count*3)],function(X) as.numeric(Get_Peptide_Scores(PEPTIDE_INPUT, X)))
      
      cat(paste(" Calculating", (length(peptable.peptide.subset)-quarter_count*3)," scores. This may take a few minutes\n"))
      TEST_Scores.filter_4 = sapply(peptable.peptide.subset[(quarter_count*3 + 1):length(peptable.peptide.subset)],function(X) as.numeric(Get_Peptide_Scores(PEPTIDE_INPUT, X)))
      
      TEST_Scores.filter = c(TEST_Scores.filter_1,TEST_Scores.filter_2,TEST_Scores.filter_3,TEST_Scores.filter_4)
      
      peptable$score = -100
      
      cat(" Sorting Scores\n")
      peptable$score[match(peptable.peptide.subset,as.character(peptable$peptide))] = TEST_Scores.filter
      
      peptable$BindLevel = as.character(peptable$BindLevel)
      peptable$BindLevel[which(peptable$BindLevel == "<= WB")] <- "Weak Binder"
      peptable$BindLevel[which(peptable$BindLevel == "<= SB")] <- "Strong Binder"
      
    }
    
    Return_TOP = 100
    
    colnames(peptable) <- c("Peptide", "Gene", "Affinity_(nM)","Bind_Level", "Peptide_Score")
    peptable$Gene <- as.character(peptable$Gene)
    peptable <- peptable[order(peptable$Peptide_Score, decreasing = T)[1: Return_TOP],]
    
    peptable <- normal(peptable)
    
    peptable <- ligandome(peptable)

    peptable$Overall_Score = round(peptable$Peptide_Score * peptable$max_norm/peptable$'Affinity_(nM)',2)
    
    Output_peptable = peptable[,c("Peptide","Gene","Peptide_Score","Overall_Score","Affinity_(nM)","max_norm","max_tissue","Normal_Ligandome","HLA")]
    
    Output_peptable$max_tissue <- as.character(Output_peptable$max_tissue) 
    Output_peptable$max_norm <- round(as.numeric(Output_peptable$max_norm),2)
    
    Output_peptable = Output_peptable[order(Output_peptable$Overall_Score,decreasing=T),]
    
    rownames(Output_peptable) <- NULL
    
    if(input_filetype == "csv")
    {
      write.csv(Output_peptable,Output_Path)
    }
    if(input_filetype == "tsv")
    {
      write.table(x = Output_peptable,file = Output_Path,quote = F,row.names = F,col.names = T,sep = "\t")
    }
    cat(paste(" Results are saved in ",out_file,"\n",sep=" "))
    cat(" Job complete\n")

} #if(CHECKS) (end)
    
    
    
    
    