#' Folder Setup
#'
#'This function creates the required folder setup in the working directory.
#'
#'@return None
#'
#'@examples Folders()
#'
#'@export

Folders <- function(){
  if (!dir.exists("0.backup")){dir.create("0.backup")
  } else {warning("0.backup exists, no need to create")}
  if(!dir.exists("1.rawreads")){dir.create("1.rawreads")
  }else{warning("1.rawreads exists, no need to create")}
  if (!dir.exists("2.mergedreads")){dir.create("2.mergedreads")
  } else {warning("2.mergedreads exists, no need to create")}
  if(!dir.exists("3.strippedreads")){dir.create("3.strippedreads")
  }else{warning("3.strippedreads, no need to create")}
  if (!dir.exists("4.pooledsamples")){dir.create("4.pooledsamples")
  } else {warning("4.pooledsamples exists, no need to create")}
  if(!dir.exists("5.OTUs")){dir.create("5.OTUs")
  }else{warning("5.OTUs exists, no need to create")}
  if (!dir.exists("6.mappings")){dir.create("6.mappings")
  } else {warning("6.mappings exists, no need to create")}
  if(!dir.exists("7.DADA2")){dir.create("7.DADA2")
  }else{warning("7.DADA2, no need to create")}
  if(!dir.exists("8.LULU")){dir.create("8.LULU")
  }else{warning("8.LULU, no need to create")}
  return("Folders set up.")
}


