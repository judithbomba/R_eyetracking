require(eyelinker)
require(dplyr)
require(plyr)
require(tidyr)
require(intervals)
require(stringr)
require(parsnip)
require(assertthat)

subjects <- c(seq(from = 27, to = 38)) #c(seq(from = 27, to = 41))
days <- c(1,2)
f <- c(1,2)

wdName <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wdName)

folderAsc <- wdName #folderAsc <- "C:/Users/judit/Desktop/Bachelor/Bachelorarbeit/BA R Project/Daten/BigSis/data_eye/ascs"
filenameRdata <- paste0(wdName,"/df_September.RData") #filenameRdata <- "C:/Users/judit/Desktop/Bachelor/Bachelorarbeit/BA R Project/Daten/BigSis/analysis/September/df_September.RData"

dryRun <- FALSE # die Namen der auszulesenden Dateien vorher einmal anzeigen lassen

load(filenameRdata) # should create df_September
behav_September <- df_September
rm(df_September)

my_dataloading <- function(fullpathFilename, filename){ #hier auch zwischen Linux und Windows trennbar (if C:/ als Pfadanfang windows, ansonsten unix/linux [hier eventuell wie in Matlab eine "isUnix" ähnliche Funktion?])
  print("loading data..")
  
  print(fullpathFilename)
  assign(filename, read.asc(fullpathFilename)) #in R: . == _ ; hätte auch read_asc schreiben können
  data <- get(filename)
  
  print(paste("loaded df ", filename))
  return(data)
}
my_downsampling <- function(data, name){ #TODO in zwei Funktionen: 1. downsampling startet mit input data, 2. data laden (generieren aus Pfaden)
  
  msg <- data$msg
  raw <- data$raw
  sac <- data$sacc
  blk <- data$blinks
  fix <- data$fix
  
  assertthat::validate_that(unique(raw$cr.info)=="...", msg = paste0("problem with cr.info in data: ", name)) #stoppt den code nicht, gibt nur error aus
  #print(paste("cr.info (sollte '...' sein): ", unique(raw$cr.info), sep = "")) # wenn nur "..." alles ok
  
  #regression
  interval <- dplyr::filter(msg, str_detect(text, "Rotation"))
  Sac <- cbind(sac$stime, sac$etime)
  
  df_pupilReg <- cbind(x = raw$xp, y = raw$yp, ps = raw$ps, t = raw$time)
  df_pupilReg <- as.data.frame(df_pupilReg)
  
  #herausfinden, wann erste Sakkade in dem Zeitraum endet (=value)
  value <- df_pupilReg %>%
    mutate(sac = t %In% Sac,
           sac.index = whichInterval(t, Sac)) %>%
    group_by(sac.index) %>%
    filter(between(t, as.numeric(interval$time[1]), as.numeric(interval$time[2]))) %>%
    filter(row_number() == n()) %>%
    .$t #TODO: Funktion, der mir Wert einer konkreten Spalte geben kann hier einfügen (tydiverse values column suchen), ich mische hier die tidy pipeline mit conventional r
  
  df_pupilReg <- df_pupilReg %>% 
    filter(between(t, value[1]+1, as.numeric(interval$time[2])))
  
  # multiple lineare Regression ps ~ x + y
  lm_model <- linear_reg() %>%
    set_engine("lm") 
  
  ps_fit <- lm_model %>%
    fit(ps~x+y, data = df_pupilReg)
  ps_fit
  b1 <- ps_fit$fit$coefficients[2] # x slope
  b2 <- ps_fit$fit$coefficients[3]
  
  print("did regression")
  # pupil size um x und y slope korrigierenmutate(ps corrected = ps - b1 * x - B2 * y)
  data_ps <- dplyr::select(raw, time, xp, yp, block, ps) %>%
    mutate(ps_corrected = ps-(b1*xp)-(b2*yp))
  
  print("corrected")
  
  #Fixation
  Fix <- cbind(fix$stime, fix$etime)
  data_ps <- data_ps %>%
    mutate(fixation = time %In% Fix,
           fixation.index = whichInterval(time, Fix))
  #Sakkaden
  Sac <- cbind(sac$stime, sac$etime)
  data_ps <- data_ps %>%
    mutate(saccade = time %In% Sac,
           saccade.index = whichInterval(time, Sac))
  #Blinzeln
  Blk <- cbind(blk$stime, blk$etime)
  data_ps <- data_ps %>%
    mutate(blink = time %In% Blk,
           blink.index = whichInterval(time, Blk))
  
  # downsampling
  data_ps_100Hz <- dplyr::filter(data_ps, row_number() %% 10 == 0)
  for (i in 1:length(data_ps_100Hz$ps)){ #TODO: rewrite as mutate! hier sonst nicht klar, was der code tut!
    #print(i)
    #print(data_ps_100Hz$ps_corrected[i])
    if (!is.na(data_ps_100Hz$ps_corrected[i])){
      if (data_ps_100Hz$ps_corrected[i] <= 0){
        data_ps_100Hz$ps_corrected[i] = 0
      } 
    }
  }
  print("sampled down")
  #messages
  # downsampling zu 100 Hz beachten
  listPatterns <- c("ITI", "calibration", "nswer", "cale", "ecognition", "word", "jitter", "entence") 
  listEventnames <- c("iti", "cal", "ans", "sca", "rec", "wrd", "jit", "sen")
  stopifnot(length(listPatterns) == length(listEventnames)) #assertion als codestop
  for (ipattern in 1:length(listPatterns)){
    eventname <- listEventnames[ipattern]
    #print(eventname)
    assign(eventname, filter(msg, str_detect(text, listPatterns[ipattern])))
    #print("assigned")
    #print(length(get(eventname)$text))
    assertthat::validate_that(length(get(eventname)$text) > 0, msg = paste0("length 0: ", name, " [skipping messages]"))
    #hier unteren Code einfügen
    nearest_time <- c() # necessary to reset this to start with clean plate
    if (length(get(eventname)$text) > 0){
      for (i in seq(from = length(get(eventname)$text), to = 1, by = -1)) {
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= get(eventname)$time[i]]) #potentiell mit summarize lösbar?
        nearest_time[i] <- data_ps_100Hz$time[x] #TODO: Alternative hierfür finden? super ineffizient!
      }
    }
    
    #get(eventname)$t_dwnspld <- nearest_time #throws error
    temp <- get(eventname)
    temp$t_dwnspld <- nearest_time
    #print(temp)
    assign(eventname, temp)
    #print(get(eventname))
    #print("ls is: ")
    #print(ls())
  }
  #print("ls out of function: ")
  #print(ls())
  #print(list(lapply(listEventnames, get)))
  msg_select <- join_all(dfs = list(iti, cal, ans, sca, rec, wrd, jit, sen), type='full') %>% #list(iti, cal, ans, sca, rec, wrd, jit, sen)
    select(t_dwnspld, text) %>%
    arrange(t_dwnspld) 
  msg_select <- dplyr::rename(msg_select, time = t_dwnspld)
  
  data_ps_100Hz <- left_join(data_ps_100Hz, msg_select, by = "time")
  
  print("messages included")
  
  #assign(paste0(name,"_ps_100Hz"), data_ps_100Hz)
  
  #print("remove data")
  #list <- c("data_ps_100Hz", paste0("asc_", name))
  #rm(list = list)
  
  #for saving data
  
  return(data_ps_100Hz)
}

#fullpathFilename <- c()
for (subject in subjects){
  for (day in days){
    for (fs in f){
      startTime <- Sys.time()
      attach <- TRUE
      
      name <- sprintf("d%ip%03if%i", day, subject, fs)
      print(name)
      
      guessedFilenameStart <- sprintf("^d%ip%03if%i.*\\.asc", day, subject, fs)
      foundFiles <- dir(path = folderAsc, pattern = guessedFilenameStart)
      
      if (length(foundFiles) != 1){
        cat(sprintf("subject: %i, day: %i === nr of files identified: %i\n",
                    subject, day, length(foundFiles)))
        for (i in 1:length(foundFiles)){ #in iFile ändern -> bei integer immer i davor, wenn nur file wäre komplexerer Inhalt
          attach <- FALSE
          print(foundFiles[i])
          name <- gsub('.{0,4}$', '', foundFiles[i])
          fullpathFilename <- paste0(folderAsc,"/",as.character(foundFiles[i]))
          filename <- paste("asc_", name, sep="") #string of filename minus last 4 char
          #print(filename)
          if (!dryRun){
            name_df <- paste0(name, "_ps_100Hz")
            filename <- paste0(name_df,".RData")
            
            data <- my_dataloading(fullpathFilename, filename) #function to load data (Windows)
            assign(name_df, my_downsampling(data, name))  #function to sample down
            
            #save data
            print("save data")
            
            save(list=name_df, file = filename)
            endTime = Sys.time()
            print(endTime-startTime)
          }
        }
        
        
      }else{
        print(paste0(folderAsc,"/",foundFiles))
        fullpathFilename <- paste0(folderAsc,"/",foundFiles)
        filename <- paste("asc_", name, sep="")
        print(filename)
        if (!dryRun){
          name_df <- paste0(name, "_ps_100Hz")
          
          data <- my_dataloading(fullpathFilename, filename) #function to load data (Windows)
          assign(name_df, my_downsampling(data, name)) #function to sample down
          
          if (day == 1 && attach ==TRUE){
            #behavioral Data anfügen
            name_ps <- paste(name, "_ps_100Hz", sep = "")
            
            #index word
            #TODO: als modulo umschreiben für word.index, nur word.contin generieren
            #mutate(word.contin%%8) dann alle 0 in der index variable = 8 setzen
            indexWord <- c()
            indexContin <- c()
            j <- 1 # 1st Word
            k <- 1 # 2nd Word
            l <- 1 # 1st ITI
            m <- 1 # 2nd ITI
            
            z <- 1 # 1st Word
            c <- 1 # ITI
            s <- 1 # 2nd Word
            ß <- 1 # ITI
            pattern1 <- as.character("(1st) New-word")
            pattern2 <- as.character("(2nd) New-word")
            pattern3 <- as.character("(1st) ITI start")
            pattern4 <- as.character("(2nd) ITI on")
            data <- get(name_ps)
            for (i in 1:length(data$text)){
              txt <- data$text[i]
              #print(txt)
              if (j == 9){j <- 1}
              if (k == 9){k <- 1}
              if (l == 9){l <- 1}
              if (m == 9){m <- 1}
              if(!is.na(txt)){
                if(grepl(pattern = pattern1, txt, fixed = TRUE)){
                  indexWord <- c(indexWord, j)
                  indexContin <- c(indexContin, z)
                  j <- j+1
                  z <- z+1
                }else if(grepl(pattern = pattern2, txt, fixed = TRUE)){
                  indexWord <- c(indexWord, k)
                  indexContin <- c(indexContin, s)
                  k <- k+1
                  s <- s+1
                }else if(grepl(pattern = pattern3, txt, fixed = TRUE)){
                  indexWord <- c(indexWord, l)
                  indexContin <- c(indexContin, c)
                  l <- l+1
                  c <- c+1
                }else if(grepl(pattern = pattern4, txt, fixed = TRUE)){
                  indexWord <- c(indexWord, m)
                  indexContin <- c(indexContin, ß)
                  m <- m+1
                  ß <- ß+1
                }else{
                  indexWord <- c(indexWord, NA)
                  indexContin <- c(indexContin, NA)
                }
              }else{
                indexWord <- c(indexWord, NA)
                indexContin <- c(indexContin, NA)
              }
            }
            data$word.index <- indexWord
            data$word.contin <- indexContin
            
            print("data indexed")
            
            # behavioral data subsetten und anfügen
            sub <- subject
            if (fs == 1){fFeedback <- "Human"}else if (fs == 2){fFeedback <- "PC"}
            behav_temp <- behav_September %>%
              subset(subject == sub) %>%
              subset(feedbackType == fFeedback) %>%
              select(congruent, correct, feedback, order, feedbackType, stimNewWord1, stimMeaningInS1, stimMeaningInS2, stimSentence1, stimSentence2) %>%
              mutate (word.contin = seq(1,88))
            
            data_join <- left_join(data, behav_temp, "word.contin")
            
            name_df <- paste(name_ps, "_behav", sep = "")
            filename <- paste(name_df, ".RData", sep = "")
            
            assign(name_df, data_join)      
            
            print(sprintf('%s assigned', name_df))
          }else{
            name_df <- sprintf("%s_ps_100Hz", name)
            filename <- paste(name_df,".RData",sep="")
          }
          
          save(list=name_df, file = filename)
          endTime = Sys.time()
          print(endTime-startTime)
          }
        }
        
        #code für Sonderfälle?
        
      }
      
  }
}
