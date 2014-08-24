description_version<-function(){
cat("\n")
cat("-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\n\n")
cat("   8888888b.  888    d8P   .d888 d8b 888     \n")
cat("   888   Y88b 888   d8P   d88P   Y8P 888     \n")
cat("   888    888 888  d8P    888        888     \n")
cat("   888   d88P 888d88K     888888 888 888888  \n")
cat("   8888888P   8888888b    888    888 888     \n")
cat("   888        888  Y88b   888    888 888     \n")
cat("   888        888   Y88b  888    888 Y88b.   \n")
cat("   888        888    Y88b 888    888   Y888  \n\n")
cat("-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\n\n")
cat(" This report was generated using PKfit v1.2.3\n")
cat(" on:-",date(),"\n")
username<-Sys.info()[['user']]
osname_version<-c(paste(Sys.info()[['sysname']],"-",Sys.info()[['version']],"\n",
                  Sys.info()[['release']],",",Sys.info()[['machine']]))
cat(" running on:",osname_version,"\n")
cat(" user id:",username,"\n\n")
cat(" PKfit is developed by Chun-ying Lee & Yung-jin Lee.\n")
cat(" contact: Yung-jin Lee <mobilepk at gmail.com> \n\n")
cat(" PKfit is under license of GPL-2|GPL-3.\n\n")
cat("..................................................")
}