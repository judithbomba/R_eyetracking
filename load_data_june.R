require(eyelinker)
require(dplyr)
require(plyr)
require(tidyr)
require(intervals)
require(stringr)
require(parsnip)

subjects <- 21 #c(seq(from = 21, to = 26))
days <- 1 #c(1,2)
f <- 2 #c(1,2)
ascFolder <- "C:/Users/judit/Desktop/Bachelor/Bachelorarbeit/BA R Project/Daten/BigSis/data_eye/ascs"

fullpathFilename <- c()

load("C:/Users/judit/Desktop/Bachelor/Bachelorarbeit/BA R Project/Daten/BigSis/analysis/June/df_June.RData")
behav_June <- df_June

for (subject in subjects){
  for (day in days){
    for (fs in f){
      startTime = Sys.time()
      name <- sprintf("d%ip%03if%i", day, subject, fs)
      print(name)
      
      guessedFilenameStart <- sprintf("^d%ip%03if%i.*\\.asc", day, subject, fs)
      foundFiles <- dir(path = ascFolder, pattern = guessedFilenameStart)
      
      if (length(foundFiles) != 1){
        cat(sprintf("subject: %i, day: %i === nr of files identified: %i\n",
                    subject, day, length(foundFiles)))
        #print(foundFiles)
        print(paste0(ascFolder,"/",foundFiles[2]))
        fullpathFilename <- c(fullpathFilename,paste0(ascFolder,"/",foundFiles[2]))
        toParse <- sprintf("asc_%s <- read.asc('%s')", name, fullpathFilename)
        print(toParse)
        eval(parse(text = toParse))
      }else{
        print(paste0(ascFolder,"/",foundFiles))
        fullpathFilename <- c(fullpathFilename,paste0(ascFolder,"/",foundFiles))
        toParse <- sprintf("asc_%s <- read.asc('%s')", name, fullpathFilename)
        print(toParse)
        eval(parse(text = toParse))
      }
      
      print("loaded df")
      
      #hier Code auslesen Variablen
      toParse <- sprintf("msg <- asc_%s$msg ", name)
      eval(parse(text = toParse))
      #msg <- asc_d1p021f1$msg
      toParse <- sprintf("raw <- asc_%s$raw ", name)
      eval(parse(text = toParse))
      #raw <- asc_d1p021f1$raw
      toParse <- sprintf("sac <- asc_%s$sacc ", name)
      eval(parse(text = toParse))
      #sac <- asc_d1p021f1$sacc
      toParse <- sprintf("blk <- asc_%s$blinks ", name)
      eval(parse(text = toParse))
      #blk <- asc_d1p021f1$blinks
      toParse <- sprintf("fix <- asc_%s$fix ", name)
      eval(parse(text = toParse))
      #fix <- asc_d1p021f1$fix
      
      print(unique(raw$cr.info)) # wenn nur "..." alles ok
      
      
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
        .$t
      
      df_pupilReg <- filter(df_pupilReg, between(t, value[1]+1, as.numeric(interval$time[2])))
      
      
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
      # mean(data_ps$ps_corrected, na.rm = T)
      # median(data_ps$ps_corrected, na.rm = T)
      # range(data_ps$ps_corrected, na.rm = T)
      
      print("corrected")
      #Fixation
      Fix <- cbind(fix$stime, fix$etime)
      data_ps <- data_ps %>%
        mutate(fixation = time %In% Fix,
               fixation.index = whichInterval(time, Fix))
      
      # downsampling
      data_ps_100Hz <- dplyr::filter(data_ps, row_number() %% 10 == 0) # von ~ 15mb auf 1.5 reduziert <3
      for (i in 1:length(data_ps_100Hz$ps)){
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
      iti <- filter(msg, str_detect(text, "ITI"))
      cal <- filter(msg, str_detect(text, "calibration"))
      ans <- filter(msg, str_detect(text, "nswer"))
      sca <- filter(msg, str_detect(text, "cale"))
      rec <- filter(msg, str_detect(text, "ecognition"))
      wrd <- filter(msg, str_detect(text, "word"))
      jit <- filter(msg, str_detect(text, "jitter"))
      sen <- filter(msg, str_detect(text, "entence"))
      
      #take the smalles nearest value 
      nearest_timeITI <- c()
      nearest_timeCAL <- c()
      nearest_timeANS <- c()
      nearest_timeSCA <- c()
      nearest_timeREC <- c()
      nearest_timeWRD <- c()
      nearest_timeJIT <- c()
      nearest_timeSEN <- c()
      
      for (i in 0:length(iti$text)) {
        #print(iti$time[i])
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= iti$time[i]])
        #print(data_ps_100Hz$time[x])
        #print(x)
        nearest_timeITI <- c(nearest_timeITI,data_ps_100Hz$time[x])
      }
      for (i in 0:length(cal$time)) {
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= cal$time[i]])
        nearest_timeCAL <- c(nearest_timeCAL,data_ps_100Hz$time[x])
      }
      for (i in 0:length(ans$time)) {
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= ans$time[i]])
        nearest_timeANS <- c(nearest_timeANS,data_ps_100Hz$time[x])
      }
      for (i in 0:length(sca$time)) {
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= sca$time[i]])
        nearest_timeSCA <- c(nearest_timeSCA,data_ps_100Hz$time[x]) 
      }
      for (i in 0:length(rec$time)) {
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= rec$time[i]])
        nearest_timeREC <- c(nearest_timeREC,data_ps_100Hz$time[x]) 
      }
      for (i in 0:length(wrd$time)) {
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= wrd$time[i]])
        nearest_timeWRD <- c(nearest_timeWRD,data_ps_100Hz$time[x]) 
      }
      for (i in 0:length(jit$time)) {
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= jit$time[i]])
        nearest_timeJIT <- c(nearest_timeJIT,data_ps_100Hz$time[x]) 
      }
      for (i in 0:length(sen$time)) {
        x <- which.max(data_ps_100Hz[data_ps_100Hz <= sen$time[i]])
        nearest_timeSEN <- c(nearest_timeSEN,data_ps_100Hz$time[x]) 
      }
      iti$t_dwnspld <- nearest_timeITI
      cal$t_dwnspld <- nearest_timeCAL
      ans$t_dwnspld <- nearest_timeANS
      sca$t_dwnspld <- nearest_timeSCA
      rec$t_dwnspld <- nearest_timeREC 
      wrd$t_dwnspld <- nearest_timeWRD 
      jit$t_dwnspld <- nearest_timeJIT
      sen$t_dwnspld <- nearest_timeSEN 
      
      msg_select <- join_all(list(iti, cal, ans, sca, rec, wrd, jit, sen), type='full') %>%
        select(t_dwnspld, text) %>%
        arrange(t_dwnspld) 
      msg_select <- dplyr::rename(msg_select, time = t_dwnspld)
      
      data_ps_100Hz <- left_join(data_ps_100Hz, msg_select, by = "time")
      
      print("messages included")
      
      
      toParse <- sprintf("%s_ps_100Hz <- data_ps_100Hz", name)
      eval(parse(text = toParse))
      
      print("remove data")
      rm(data_ps_100Hz)
      toParse <- sprintf("rm(asc_%s)", name)
      eval(parse(text = toParse))
      
      # 
      # behavioral data anfügen
      #
      if (day == 1){
        name_ps <- sprintf("%s_ps_100Hz", name)
        
        #index word
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
        behav_temp <- behav_June %>%
          subset(subject == sub) %>%
          subset(feedbackType == fFeedback) %>%
          select(congruent, correct, feedback, order, feedbackType, stimNewWord1, stimMeaningInS1, stimMeaningInS2, stimSentence1, stimSentence2) %>%
          mutate (word.contin = seq(1,88))
        
        data_join <- left_join(data, behav_temp, "word.contin")
        
        name_df <- sprintf('%s_behav',name_ps)
        filename <- paste(name_df,".RData",sep="")
        
        assign(name_df, data_join)      
        
        print(sprintf('%s assigned',name_df))
      }
      
      save(list=name_df, file = filename)
      endTime = Sys.time()
      print(endTime-startTime)
      
    }
  }
}
