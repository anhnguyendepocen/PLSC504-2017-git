#################################################
# Presentation order for PLSC 504 - December 2017
#
# Six presentations on Tuesday (12/5), seven on 
# Thursday (12/7).
#
# File created 12/1/2017
#################################################

# Class members (alphabetically):

peeps <- c("Brandon Bolte","Jinhyuk Jang","Boyoon Lee",
           "Steven Morgan","Rosemary Pang","Joseph Phillips",
           "Connor Somgynari","Michael Soules","Kim Tran",
           "Aubrey Waddick","Mikaela Westhoff","Samuel Wilmer",
           "Omer Yalcin")

# Data:

df <- data.frame(Name = peeps)

# Flag Tuesday volunteers (Lee, Phillips, Waddick, Westhoff):

df$Tuesdummy <- c(0,0,1,0,0,1,0,0,0,1,1,0,0)

# Draw pseudo-random numbers; the non-Tuesday class members
# with the highest two values are assigned to Tuesday:

set.seed(7222009)
df$rando <- rnorm(nrow(df),0,1)
df <- df[order(df$Tuesdummy, df$rando),]
df$Tuesdummy[8:13] <- 1 

# Now re-sort within Tuesday/Thursday groups to determine
# presentation order:

df <- df[order(-df$Tuesdummy, df$rando),]
df$Day <- ifelse(df$Tuesdummy==1,paste("Tuesday"),paste("Thursday"))
df$rando <- NULL
df$Tuesdummy <- NULL

df

# fin