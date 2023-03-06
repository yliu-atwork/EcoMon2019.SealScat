# Functions needed for calculation of Multifrequency Indicator and Areal Backscatter Coefficient from subset Sv data
# and for comparison of or ABC and eDNA vertebrate band intensity values
# Developer: Skylar Gering 08/2021

library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)

calc_MFI <- function(ctd_depth, bad_fq = character(), delta = 40, global_norm = FALSE){
  # Calculate MFI from subsetted Sv data in JSON file (from segments_subsets_mfi.py) and classifying functions
  # Inputs: ctd_depth - list of subsets (arrays of Sv data), one for each frequency at a given depth; bad_fq - list
  # of string frequencies (i.e. "18000") that should not be used; delta - integer from MFI formula (see Trenkel, Verena ,and Laurent Berger.);
  # global_norm - boolean determines local (one frequency) or global (all frequencies) normalization of data
  # Outputs: one array of MFI values (same dimensions as each of the frequency Sv arrays)
  
  ctd_depth <- ctd_depth[names(ctd_depth) %in% bad_fq == FALSE]
  # create 3D array of Sv data - one slice per frequency
  nd = ncol(ctd_depth[[1]]) # number of depth bins
  np = nrow(ctd_depth[[1]]) # number of pings
  nf = length(ctd_depth) # number of frequencies
  
  if (nf > 2){
    Sv_f = array(0, dim=(c(np, nd, nf)))
    
    for (i in 1:nf){
      Sv_f[, , i] = ctd_depth[[i]]
    }
    # linear data
    sv_f = 10^(Sv_f/10)
    
    # frequencies
    f <- unlist(lapply(names(ctd_depth), function(x) strtoi(x)/1000))
    e_f = 1/f
    
    # all unique combinations of indices
    f_idx = combinations(n=nf, r=2, set=F, repeats.allowed=F)
    d_f = cbind(f_idx, 0)
    
    # Distance function
    for (i in 1:length(f_idx[,1])) {
      d_f[i,3] = 1-exp(-abs(f[f_idx[i,1]]-f[f_idx[i,2]])/delta)
    }
    
    #scaled linear data
    D_f = sv_f
    
    if (!global_norm) { # scale to local (per frequency) maximum
      for (fq in 1:nf) {
        D_f[,,fq] = sv_f[,,fq]/max(sv_f[,,fq], na.rm=TRUE)
      }
    }
    if (global_norm) { # scale to global maximum
      D_f = sv_f/max(sv_f, na.rm=TRUE)
    }
    
    # calculate the MFI values
    MFI = array(0, dim=(c(np, nd)))
    num = array(0, dim=(c(np, nd)))
    den = array(0, dim=(c(np, nd)))
    for (i in 1:(nf-1)) {
      for (l in (i+1):nf) {
        num = num+
          d_f[which(d_f[,1] == i & d_f[,2] == l),3]*D_f[,,i]*D_f[,,l]*e_f[i]*e_f[l]
        den = den+D_f[,,i]*D_f[,,l]*e_f[i]*e_f[l]
      }
    }
    return(((num/den)-0.4)/0.6)
  }
  else {
    print("MFI cannot be calculated with less than 3 frequencies.")
    return (NULL)
  }
}

classify_MFI <- function(MFI_arr){
  # Separate MFI values into 4 catagories as defined by Trenkel, Verena M., and Laurent Berger paper
  # "A fisheries acoustic multi-frequency indicator to inform on large scale spatial patterns of aquatic pelagic ecosystems."
  # Inputs: MFI_arr - array of MFI data
  # Outputs: array with values replaced with catagory name
  
  bins <- c(-Inf, 0.4, 0.6, 0.8, Inf)  # some values are very slightly smaller that 0
  catagories <- c("Swimbladder Fish", "Small Bubbles", "Zooplankton", "Non-swimbladder Fish")
  class_MFI <- cut(MFI_arr, breaks = bins, labels = catagories)
  return(class_MFI)
}

isclose <- function(a, b, abs_tol){
  # Returns True if a and b are within abs_tol of each other, False otherwise
  return(abs(a-b) <= abs_tol)
}

abc_df <- function(data, cast)
  # takes in ABC JSON file from mfi_masks_abc.py python code and cast number and unpacks data 
  # from that cast and formats into data frame
  # Inputs: data - abc JSON data, cast - integer cast number
  # Outputs: data frame with abc data and cast/depth
  # cat("cast = ", cast)
  df <- do.call(rbind, data[[cast]]) # selects one cast's data
  # df <- do.call(rbind, data[[names(data)]])
  if(!is.null(dim(df))){ # if there is data for given cast
    df <- cbind(rownames(df), data.frame(df, row.names=NULL)) %>% mutate(Cast = cast) # makes row names a column
    # problem is probably in previous line maybe in the mutate?
    print(df)
    colnames(df) <- list("Depth", "ABC", "Cast")
  }
  return(df)


approx_depth_mutate <- function(edna, df){
  # For one cast add eDNA.Depth column to ABC data frame by matching eDNA depth and ABC depth - uses isclose function as depths
  # recorded by echosounder and CTD rosette can be slightly different 
  # Inputs: edna - data frame from reading in eDNA data with the depth column called "eDNA.Depth", abc data frame from reading in 
  # with abc_df with deoth column called "Depth"
  # Ouput: updated abc data frame for one cast with added column called "eDNA.Depth"
  return (mutate(df, eDNA.Depth = unlist(lapply(df$Depth, function(x) edna[isclose(strtoi(x), edna$eDNA.Depth, 5), ]$eDNA.Depth))))
}

mask_MFI <- function(MFI_arr, Sv_arr, ranges){
  # Creates a mask that preserves data within the given ranges for the MFI_arr and applies that mask to the Sv_arr
  # Inputs: MFI_arr - array of MFI data;  Sv_arr - array of Sv data; ranges - list of tuples where the first number is the
  # bottom of the range and the second is the top
  # Output: Sv_array with numbers outside of ranges set to -999
  # Note:  MFI_arr and Sv_arr need to have the same dimensions
  mask <- Reduce("|", lapply(ranges, function(r) between(MFI_arr, r[[1]], r[[2]])))
  Sv_masked <- as.double(mask) * Sv_arr
  Sv_masked <- replace(Sv_masked, Sv_masked == 0, -999)
  return(Sv_masked)
}

calc_ABC <- function(masked_Sv, depth_bounds){
  # takes a masked_Sv array and calculates the Areal Backscattering Coefficient (ABC)
  # Inputs: masked_Sv - Sv array, depth_bounds - 1D array with depth ticks for the Sv array
  # Output: float with ABC value
  mean_val <- sum(10^(masked_Sv/10))/length(masked_Sv)
  bin_thickness <- (max(depth_bounds) - min(depth_bounds))/length(depth_bounds)
  return(mean_val*bin_thickness)
}

abc_edna_violin <- function(df){
  # plotting function for ABC value vs the eDNA band intensity rank violin plot
  # Need a data frame with both ABC and X.12S.final columns
  p <- ggplot(df, aes(X.12S.final, log(ABC), fill = factor(X.12S.rank))) + 
    geom_violin(width=0.8) + geom_point() + 
    geom_boxplot(width=0.2, color="black", alpha=0.2) + ggtitle("Fish ABC vs 12S PCR eDNA Rank") +
    xlab("Vertebrate Band Intensity") + ylab("log(ABC) for Fish Scattering") + theme_bw() + 
    theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = 0.5), 
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))
  
  return (p)
}

abc_edna_scatter <- function(df){
  # plotting function for ABC value vs the eDNA band intensity rank scatter plot
  # Need a data frame with both ABC and X.12S.final columns
  p <- ggplot(df, aes(X.12S.rank, log_ABC, color = factor(X.12S.rank))) +
    geom_point() + theme_bw() + 
    theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = 0.5), 
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16)) + 
    ggtitle("Fish ABC vs 12S PCR eDNA Rank") +
    xlab("12S PCR eDNA Rank") + ylab("log of ABC for Fish Scattering") + 
    scale_x_continuous(labels= c("No", "EL", "L", "M", "H"))
  
  return (p)
}

read_eDNA <- function(edna_data, transducer_depth = 5, filtration_vol){
  # Reads in very specific eDNA CSV - will need to be edited to match other formats
  # Select specific columns from eDNA CSV and create new column called X.12s.rank which gives a number rank from 0-4 for eDNA 12s band intensity
  edna_data %>% 
    select(Station, Cast, Lat, Long, Filtration.Volume, Sampling.Depth.Meter, Sampling.Depth.Type, DNA.Conc, X.12S.final) %>%
    mutate(X.12S.rank = case_when(
      endsWith(X.12S.final, "No") ~ 0,
      endsWith(X.12S.final, "EL") ~ 1,
      endsWith(X.12S.final, "L") ~ 2,
      endsWith(X.12S.final, "M") ~ 3, 
      endsWith(X.12S.final, "H") ~4)) %>% na.omit() %>%
    transform(X.12S.final=factor(X.12S.final,levels=c("No", "EL", "L", "M", "H"))) %>%
    subset(Sampling.Depth.Meter > transducer_depth & Filtration.Volume == filtration_vol) -> edna_df
  
  return(edna_df)
}
  

box_plot <- function(df){
  # plotting function for ABC value vs the eDNA band intensity rank scatter plot
  # Need a data frame with both ABC and X.12S.final columns
  p <- ggplot(df, aes(X.12S.rank, log_ABC, color = factor(X.12S.rank), fill = factor(X.12S.rank))) +
    geom_boxplot(width=0.2, color="black", alpha=0.2) + ggtitle("Fish ABC vs 12S PCR eDNA Rank") +
    #geom_point() + 
    xlab("Vertebrate Band Intensity") + ylab("log(ABC) for Fish Scattering") + theme_bw() + 
    geom_jitter(shape=16, width = 0.1, height = 0.1) + geom_text(label=df$Cast_Depth, check_overlap = TRUE, size=4, nudge_x = 0.15) +
    theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = 0.5), 
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16)) +
    scale_x_continuous(labels= c("No", "EL", "L", "M", "H"))
  
  return (p)
}

scatter_plot <- function(df, x, y){
  p <- ggplot(df, aes({{x}}, {{y}})) +
    geom_point() +
    geom_smooth(method = lm) +
    stat_cor(method = "pearson")
  return(p)
}
