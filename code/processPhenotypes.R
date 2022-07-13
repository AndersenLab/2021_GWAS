require(tidyverse)
require(easyXpress)
require(RColorBrewer)
require(cowplot)
require(magrittr)
require(ggbeeswarm)
require(ggrepel)
require(grid)
require(gridExtra)
require(sommer)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))
today <- format(Sys.Date(), "%Y%m%d")
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = na.rm, 
                         ...)
  H <- 1.5 * stats::IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
cleanData <- function(raw.data){
  
  # Remove Low n and Censored Wells
  n.drug.df <- raw.data %>%
    dplyr::filter(n < 30,
                  n > 5,
                  is.na(well_censor))

  
  # Outlier Removal
  outlier.removed.drug.df <- n.drug.df %>%
    dplyr::group_by(drug, strain) %>%
    dplyr::mutate(median_wormlength_um = remove_outliers(median_wormlength_um)) %>%
    dplyr::filter(!is.na(median_wormlength_um))
  
  # Control Delta
  control_values <- outlier.removed.drug.df %>%
    dplyr::filter(drug %in% c("DMSO", "Water")) %>% # filter to control wells
    dplyr::group_by(strain, Assay_Bleach) %>%
    dplyr::mutate(control_pheno = mean(median_wormlength_um)) %>% # get mean control value for each trait and strain
    dplyr::distinct(control_pheno, strain, Assay_Bleach) # make it neat
  
  delta <- outlier.removed.drug.df %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., control_values) %>% # join control values for each trait
    dplyr::filter(!drug %in% c("DMSO", "Water")) %>% # filter on just exposed wells
    dplyr::mutate(median_wormlength_um_delta = median_wormlength_um - control_pheno)
  
  
  # Regression
  if(nrow(delta) > 0){
    median.reg <- lm(data = delta, formula = median_wormlength_um_delta ~ Assay_Bleach)
    delta <- delta[as.numeric(names(median.reg$residuals)),]
    
    coefs <- lm(data = delta, formula = median_wormlength_um_delta ~ Assay_Bleach - 1)
    
    bleach.coefs <- data.frame(gsub(names(coefs$coefficients), pattern = "Assay_Bleach", replacement = ""), 
               as.numeric(coefs$coefficients))
    colnames(bleach.coefs) <- c("Assay_Bleach","AssayBleach_Effect")
    regressed <- cbind(delta, median.reg$residuals) %>%
      dplyr::rename(median_wormlength_um_delta_reg = `median.reg$residuals`) %>%
      dplyr::full_join(., bleach.coefs)
  } else {
    regressed <- delta
  }
  return(regressed)
}

tx.classes <- data.table::fread(file = "data/tx.classes.csv") %>%
  dplyr::rename(trait = drug)
processData <- function(data, plot = FALSE){
  
  cleaned.data <- data %>%
    cleanData()
  
  # cleaned.data.plot <- data.cleaning.plot(cleaned.data[[2]], "Data Cleaning")
  
  if(nrow(cleaned.data) == 0){
    print("All data filtered from well.n or outliers!")
  } else if (plot == FALSE){
    print("Skipping plotting")
    return(cleaned.data)
  } else{
    strain.pal <- rep("black", length(unique(cleaned.data$strain)))
    names(strain.pal) <- unique(cleaned.data$strain)
    strain.pal[names(strain.pal) %in% "PD1074"] <- "orange"
    strain.pal[names(strain.pal) %in% "CB4856"] <- "blue"
    strain.pal[names(strain.pal) %in% "JU775"] <- "green"
    strain.pal[names(strain.pal) %in% "MY16"] <- "#67697C"
    strain.pal[names(strain.pal) %in% "ECA36"] <- "#5A0C13"
    strain.pal[names(strain.pal) %in% "XZ1516"] <- "purple"
    strain.pal[names(strain.pal) %in% "RC301"] <- "#627264"
    strain.pal[names(strain.pal) %in% "ECA396"] <- "#C51B29"
    strain.pal[names(strain.pal) %in% "CB4855"] <- "#a37000"
    strain.pal[names(strain.pal) %in% "DL238"] <- "red"
    
    
    
    
    strain.distribution.plot <- cleaned.data %>%
      dplyr::group_by(strain) %>%
      dplyr::summarise(norm.length = mean(median_wormlength_um_delta_reg),
                       sd.norm.length = sd(median_wormlength_um_delta_reg)) %>%
      dplyr::mutate(label = if_else(strain %in% c("PD1074","CB4856","JU775","MY16","ECA36","XZ1516","RC301","ECA396","CB4855","DL238"), 
                                    true = strain, 
                                    false = ""))
    
    strain.distribution.plot.2 <- ggplot(strain.distribution.plot) +
      theme_bw(base_size = 8) +
      geom_pointrange(mapping = aes(x = reorder(strain,norm.length), 
                                    y = norm.length,
                                    ymin = norm.length - sd.norm.length,
                                    ymax = norm.length + sd.norm.length,
                                    colour = strain)) +
      theme(axis.text.x = element_blank(),
            legend.position = "none",
            panel.grid.minor = element_blank()) + 
      geom_label_repel(aes(x = reorder(strain,norm.length),
                          y = norm.length,
                          label = label), size = 2, max.overlaps = Inf) + 
      scale_colour_manual(values = strain.pal, name = "Assay") + 
      labs(x = "Strain", 
           y = "Bleach-Regressed Median Worm Length (um)", title = unique(cleaned.data$drug))
  

    ggsave(strain.distribution.plot.2, filename = paste0("plots/",today,"_",unique(cleaned.data$drug),"_reg_combined_GWA_status.png"),
           width = 10, height = 5)
    
    
    
    
    
    strain.distribution.plot <- cleaned.data %>%
      ggplot(.) +
      theme_bw(base_size = 9) +
      geom_point(mapping = aes(x = reorder(strain,median_wormlength_um), 
                               y = median_wormlength_um, 
                               shape = as.factor(bleach),
                               colour = Metadata_Experiment),
                 position = position_quasirandom(width = 0.1)) +
      geom_boxplot(mapping = aes(x = reorder(strain,median_wormlength_um), 
                                 y = median_wormlength_um,
                                 fill = strain), 
                   alpha = 0.7) +
      scale_fill_manual(values = strain.pal, guide = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "bottom",
            panel.grid.minor = element_blank()) +
      scale_colour_brewer(palette = "Set3", name = "Assay") +
      scale_shape_manual(values = c(15,16,17), name = "Bleach") +
      labs(x = "Strain", y = "Median Worm Length (um)", title = unique(cleaned.data$drug))
    
    ggsave(strain.distribution.plot, filename = paste0("plots/",today,"_",unique(cleaned.data$drug),"_combined_GWA_status.png"),
           width = 10, height = 4)
    
    return(cleaned.data)
    
  }
  
}
gwa.metadata <- data.table::fread("data/gwa.doses.csv")


#################
#### GWA 01B ####
#################

GWA01B.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20210902_GWA01B/",
                                        design = F,
                                        rdafile = "CellProfiler-Analysis_20210902_GWA01B_data_1631050578SJW.RData") %>%
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA01B.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20210902_GWA01B/design/20210902_GWA01B_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA01B")

processed.GWA01B.CP.dat <- GWA01B.CP.dat %>%
  dplyr::full_join(., GWA01B.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)


#################
#### GWA 02C ####
#################

GWA02C.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20210826_GWA02C/",
                                        design = F,
                                        rdafile = "CellProfiler-Analysis_20210826_GWA02C_data_1631036380SJW.RData") %>%
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA02C.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20210826_GWA02C/design/20210826_GWA02C_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA02C")

processed.GWA02C.CP.dat <- GWA02C.CP.dat %>%
  dplyr::full_join(., GWA02C.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)


################
#### GWA 04 ####
################

GWA04.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20211008_GWA04/",
                                        design = F,
                                        rdafile = "CellProfiler-Analysis_20211008_GWA04_data_1633724904SJW.RData") %>%
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA04.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20211008_GWA04/design/20211008_GWA04_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA04")

processed.GWA04.CP.dat <- GWA04.CP.dat %>%
  dplyr::full_join(., GWA04.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)



################
#### GWA 05 ####
################

GWA05.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20211022_GWA05/",
                                       design = F,
                                       rdafile = "CellProfiler-Analysis_20211022_GWA05_data_1635029338_JWoct27_2.RData") %>% # "CellProfiler-Analysis_20211022_GWA05_data_1635029338JWoct23.RData" | "cell_profiler_numeric.RData"
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA05.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20211022_GWA05/design/20211022_GWA05_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA05")

processed.GWA05.CP.dat <- GWA05.CP.dat %>%
  dplyr::full_join(., GWA05.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)


################
#### GWA 06 ####
################

GWA06.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20211028_GWA06/",
                                       design = F,
                                       rdafile = "CellProfiler-Analysis_20211028_GWA06_data_1635619989_JWoct30.RData") %>% # "CellProfiler-Analysis_20211022_GWA05_data_1635029338JWoct23.RData" | "cell_profiler_numeric.RData"
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA06.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20211028_GWA06/design/20211028_GWA06_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA06")

processed.GWA06.CP.dat <- GWA06.CP.dat %>%
  dplyr::full_join(., GWA06.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)

################
#### GWA 07 ####
################

GWA07.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20211112_GWA07/",
                                       design = F,
                                       rdafile = "CellProfiler-Analysis_20211112_GWA07_data_1636987645SJW.RData") %>%
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA07.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20211112_GWA07/design/20211112_GWA07_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA07")

processed.GWA07.CP.dat <- GWA07.CP.dat %>%
  dplyr::full_join(., GWA07.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)


################
#### GWA 08 ####
################

GWA08.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20211118_GWA08/",
                                       design = F,
                                       rdafile = "CellProfiler-Analysis_20211118_GWA08_data_1637468337SJW.RData") %>%
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA08.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20211118_GWA08/design/20211118_GWA08_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA08")

processed.GWA08.CP.dat <- GWA08.CP.dat %>%
  dplyr::full_join(., GWA08.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)


#################
#### GWA 03B ####
#################

GWA03B.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20211210_GWA03B/",
                                       design = F,
                                       rdafile = "CellProfiler-Analysis_20211210_GWA03B_data_1639241561SJW.RData") %>%
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA03B.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20211210_GWA03B/design/20211210_GWA03B_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA03B")

processed.GWA03B.CP.dat <- GWA03B.CP.dat %>%
  dplyr::full_join(., GWA03B.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)


################
#### GWA 09 ####
################

GWA09.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20220128_GWA09/",
                                        design = F,
                                        rdafile = "CellProfiler-Analysis_20220128_GWA09_data_1643740893SJW.RData") %>%
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA09.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20220128_GWA09/design/20220201_GWA09_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA09")

processed.GWA09.CP.dat <- GWA09.CP.dat %>%
  dplyr::full_join(., GWA09.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)


################
#### GWA 10 ####
################

GWA10.CP.dat <- easyXpress::readXpress(filedir = "~/Dropbox/HTA_ImageXpress/Projects/20220203_GWA10/",
                                       design = F,
                                       rdafile = "CellProfiler-Analysis_20220203_GWA10_data_1644290191SJW.RData") %>%
  dplyr::filter(Worm_Length > 50) %>%
  easyXpress::modelSelection(.) %>%
  easyXpress::edgeFlag(.) %>%
  easyXpress::setFlags(.)

GWA10.design <- data.table::fread("~/Dropbox/HTA_ImageXpress/Projects/20220203_GWA10/design/20220203_GWA10_design.csv") %>%
  dplyr::mutate(Metadata_Experiment = "GWA10")

processed.GWA10.CP.dat <- GWA10.CP.dat %>%
  dplyr::full_join(., GWA10.design) %>%
  easyXpress::process(., Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain)



########################
#### COMBINE ASSAYS ####
########################



## combining all assays ##
prelim.GWA <- processed.GWA01B.CP.dat$summarized_processed %>%
  dplyr::full_join(., processed.GWA02C.CP.dat$summarized_processed) %>%
  dplyr::full_join(., processed.GWA04.CP.dat$summarized_processed) %>%
  dplyr::full_join(., processed.GWA05.CP.dat$summarized_processed) %>%
  dplyr::full_join(., processed.GWA06.CP.dat$summarized_processed) %>%
  dplyr::full_join(., processed.GWA07.CP.dat$summarized_processed) %>%
  dplyr::full_join(., processed.GWA08.CP.dat$summarized_processed) %>%
  dplyr::full_join(., processed.GWA03B.CP.dat$summarized_processed) %>%
  dplyr::full_join(., processed.GWA09.CP.dat$summarized_processed) %>%
  dplyr::full_join(., processed.GWA10.CP.dat$summarized_processed) %>%
  tidyr::unite("Assay_Bleach", Metadata_Experiment:bleach, remove = F)


## calculating well coefficients of variation of well census ##
well.stats <- prelim.GWA %>%
  dplyr::select(strain, drug, Metadata_Well, Assay_Bleach, bleach, n) %>%
  # dplyr::filter(!is.na(n)) %>%
  dplyr::group_by(Assay_Bleach, strain) %>%
  dplyr::summarise(num.wells = n(),
                   perc.wells = (num.wells/232)*100, # 4 wells * 58 compounds = 232
                   mean.n = mean(n),
                   cv.n = var(n)/mean(n))

## filtering out wells with CV.n greater than 5 ##
prelim.GWA.filtered <- prelim.GWA %>%
  dplyr::left_join(., well.stats) %>%
  dplyr::mutate(cv.flag = if_else(cv.n > 5, true = "CV Flag", false = "No CV Flag")) %>%
  dplyr::filter(!is.na(drug),
                cv.flag == "No CV Flag") 

grouped.prelim.GWA <- prelim.GWA.filtered %>%
  dplyr::group_by(drug) %>%
  tidyr::nest()


## ROUND 1: processing data to identify assays with high block effects ##
# Process Data and Plot Strain Distributions
assignControl <- function(grouped.data, grouped.drug){
  drug.specific.control <- unique(grouped.data$diluent)
  renamed.drug.data <- grouped.data %>%
    dplyr::mutate(drug = grouped.drug)
  if(drug.specific.control == "DMSO"){
    control.included <- prelim.GWA.filtered %>%
      dplyr::filter(drug == "DMSO") %>%
      dplyr::full_join(renamed.drug.data)
  } else if (drug.specific.control == "Water"){
    control.included <- prelim.GWA.filtered %>%
      dplyr::filter(drug == "Water") %>%
      dplyr::full_join(renamed.drug.data)
  } else {
    print("No control data; retry")
    next
  }
  return(control.included)
}
control.assigned.GWA <- purrr::map2(grouped.prelim.GWA$data,
                                    grouped.prelim.GWA$drug,
                                    assignControl)

all.processed.GWA.data <- purrr::map2(control.assigned.GWA, 
                                      FALSE,
                                      processData)
GWA.df <- all.processed.GWA.data %>% 
  purrr::keep(function(x) is.data.frame(x)) %>%
  Reduce(rbind, .)

## calculating block effect coefficients
assay_effect_coefs <- GWA.df %>%
  dplyr::select(Metadata_Experiment, bleach, drug, Assay_Bleach, AssayBleach_Effect) %>%
  dplyr::full_join(., gwa.metadata) %>%
  dplyr::filter(!drug %in% c("Water","DMSO")) %>%
  dplyr::distinct()

included.assay.bleaches <- assay_effect_coefs %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(mean.block.effect  = mean(AssayBleach_Effect), 
                   sd.block.effect = sd(AssayBleach_Effect),
                   upper = mean.block.effect + 1.5*(sd.block.effect),
                   lower = mean.block.effect - 1.5*(sd.block.effect)) %>%
  dplyr::select(-mean.block.effect, -sd.block.effect) %>%
  dplyr::left_join(.,assay_effect_coefs) %>%
  dplyr::filter(AssayBleach_Effect > lower,
                AssayBleach_Effect < upper)


ggplot(included.assay.bleaches, mapping = aes(x = Metadata_Experiment, y = AssayBleach_Effect)) + 
  theme_bw() + 
  geom_jitter(width = 0.2, height = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(.~drug) + 
  # scale_color_manual(values = c("darkred","black")) +
  labs(y = "AssayBleach Regression Coefficient",
       x = "Assay")




## ROUND 2: processing data to identify assays with high block effects ##

## filtering out wells with CV.n greater than 5 and separating out control data ##
rd2.control.data <- prelim.GWA %>%
  dplyr::left_join(., well.stats) %>%
  dplyr::mutate(cv.flag = if_else(cv.n > 5, true = "CV Flag", false = "No CV Flag")) %>%
  dplyr::filter(!is.na(drug),
                cv.flag == "No CV Flag") %>%
  dplyr::filter(drug %in% c("Water","DMSO"))

## filtering out wells with CV.n greater than 5, ##
## separating out exposure data, ##
## and filtering to low block effect assays/bleaches ##

prelim.GWA.filtered.2 <- prelim.GWA %>%
  dplyr::left_join(., well.stats) %>%
  dplyr::mutate(cv.flag = if_else(cv.n > 5, true = "CV Flag", false = "No CV Flag")) %>%
  dplyr::filter(!is.na(drug),
                cv.flag == "No CV Flag") %>%
  dplyr::filter(!drug %in% c("Water","DMSO")) %>%
  dplyr::right_join(.,included.assay.bleaches %>%
                      dplyr::select(drug, Assay_Bleach)) %>%
  dplyr::bind_rows(rd2.control.data)

grouped.prelim.GWA.2 <- prelim.GWA.filtered.2 %>%
  dplyr::group_by(drug) %>%
  tidyr::nest()

assignControl.postblockeffects <- function(grouped.data, grouped.drug){
  drug.specific.control <- unique(grouped.data$diluent)
  renamed.drug.data <- grouped.data %>%
    dplyr::mutate(drug = grouped.drug)
  if(drug.specific.control == "DMSO"){
    control.included <- prelim.GWA.filtered.2 %>%
      dplyr::filter(drug == "DMSO") %>%
      dplyr::full_join(renamed.drug.data)
  } else if (drug.specific.control == "Water"){
    control.included <- prelim.GWA.filtered.2 %>%
      dplyr::filter(drug == "Water") %>%
      dplyr::full_join(renamed.drug.data)
  } else {
    print("No control data; retry")
    next
  }
  return(control.included)
}
control.assigned.GWA.2 <- purrr::map2(grouped.prelim.GWA.2$data,
                                      grouped.prelim.GWA.2$drug,
                                      assignControl)

all.processed.GWA.data.2 <- purrr::map2(control.assigned.GWA.2,
                                      FALSE,
                                      processData)
GWA.df.2 <- all.processed.GWA.data.2 %>% 
  purrr::keep(function(x) is.data.frame(x)) %>%
  Reduce(rbind, .)


# Create trait files for each exposure for NemaScan
pre.traitfiles <- GWA.df.2 %>%
  dplyr::select(strain, drug, median_wormlength_um_delta_reg) %>%
  dplyr::group_by(strain, drug) %>%
  dplyr::summarise(mean = mean(median_wormlength_um_delta_reg)) %>%
  dplyr::mutate(drug = gsub(drug, pattern = " ", replacement = "_")) %>%
  dplyr::mutate(drug = gsub(drug, pattern = "-", replacement = "_")) %>%
  dplyr::mutate(strain = if_else(strain == "PD1074", true = "N2", false = strain)) %>%
  dplyr::group_by(drug) %>%
  tidyr::nest()

# Output all trait files
purrr::map2(.x = pre.traitfiles$data,
            .y = pre.traitfiles$drug,
            .f = function(trait.data, trait){
              colnames(trait.data) <- c("strain",trait)
              write_tsv(trait.data, file = paste0("data/",trait,"_traitfile.tsv"))
              }
            )

unique(GWA.df.2$drug) %in% unique(gwa.metadata$drug)

toxicant.traitfile <- GWA.df.2 %>%
  dplyr::select(strain, drug, median_wormlength_um_delta_reg) %>%
  dplyr::group_by(strain, drug) %>%
  dplyr::summarise(mean = mean(median_wormlength_um_delta_reg)) %>%
  dplyr::left_join(., gwa.metadata) %>%
  dplyr::filter(class == "Toxicant",
                drug != "Bacterial dilution") %>%
  dplyr::select(strain, drug, mean) %>%
  tidyr::pivot_wider(names_from = drug, values_from = mean) %>%
  dplyr::mutate(strain = if_else(strain == "PD1074", true = "N2", false = strain))
write_tsv(toxicant.traitfile, file = "output/toxicant.tratifile.tsv")

anthelmintic.traitfile <- GWA.df.2 %>%
  dplyr::select(strain, drug, median_wormlength_um_delta_reg) %>%
  dplyr::group_by(strain, drug) %>%
  dplyr::summarise(mean = mean(median_wormlength_um_delta_reg)) %>%
  dplyr::left_join(., gwa.metadata) %>%
  dplyr::filter(class == "Anthelmintic") %>%
  dplyr::select(strain, drug, mean) %>%
  tidyr::pivot_wider(names_from = drug, values_from = mean)
write_tsv(anthelmintic.traitfile, file = "output/anthelmintic.traitfile.tsv")


assay_raw_wormlength_dists <- GWA.df.2 %>%
  dplyr::select(strain, Metadata_Experiment, bleach, drug, median_wormlength_um) %>%
  # dplyr::distinct() %>%
  ggplot(., mapping = aes(x = Metadata_Experiment, y = median_wormlength_um)) + 
  theme_bw() + 
  geom_violin() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(.~drug, scales = "free_y")
ggsave(assay_raw_wormlength_dists, filename = paste0("plots/assay_raw_wormlength_dists.", today, ".png"), width = 15, height = 8)

assay_reg_wormlength_dists <- GWA.df %>%
  dplyr::select(Metadata_Experiment, bleach, drug, median_wormlength_um_delta_reg) %>%
  # dplyr::distinct() %>%
  ggplot(., mapping = aes(x = Metadata_Experiment, y = median_wormlength_um_delta_reg)) + 
  theme_bw() + 
  geom_violin() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(.~drug, scales = "free_y")
ggsave(assay_reg_wormlength_dists, filename = paste0("plots/assay_reg_wormlength_dists.", today, ".png"), width = 15, height = 8)


# Phenotype distributions 
GWA.df %>%
  dplyr::full_join(., gwa.metadata) %>%
  dplyr::filter(class == "Anthelmintic") %>%
  dplyr::select(strain, drug, median_wormlength_um_delta_reg) %>%
  dplyr::group_by(strain, drug) %>%
  dplyr::summarise(mean = mean(median_wormlength_um_delta_reg)) %>%
  dplyr::mutate(drug = if_else(drug == "Nickel dichloride", true = "Nickel chloride", false = drug),
                drug = if_else(drug == "Cadmium dichloride", true = "Cadmium chloride", false = drug),
                drug = if_else(drug == "Zinc dichloride", true = "Zinc chloride", false = drug)) %>%
  ggplot(., mapping = aes(x = mean)) + 
  theme_bw() + 
  geom_density(adjust = 0.5) + 
  facet_wrap(.~drug) + 
  theme(panel.grid = element_blank())


# Phenotype distributions filled by assay
GWA.df %>%
  dplyr::mutate(drug = if_else(drug == "Nickel dichloride", true = "Nickel chloride", false = drug),
                drug = if_else(drug == "Cadmium dichloride", true = "Cadmium chloride", false = drug),
                drug = if_else(drug == "Zinc dichloride", true = "Zinc chloride", false = drug)) %>%
  ggplot(., mapping = aes(x = median_wormlength_um_delta_reg, fill = Metadata_Experiment)) + 
  theme_bw() +
  facet_wrap(drug~., scales = "free") + 
  geom_density(adjust = 0.75, alpha = 0.1) + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_blank())


toxin.sentinel.df <- GWA.df %>%
  dplyr::filter(strain %in% c("PD1074","CB4856","MY16","JU775")) %>%
  dplyr::left_join(., gwa.metadata)
strain.pal <- rep("black", length(unique(toxin.sentinel.df$strain)))
names(strain.pal) <- unique(toxin.sentinel.df$strain)
strain.pal[names(strain.pal) %in% "PD1074"] <- "orange"
strain.pal[names(strain.pal) %in% "CB4856"] <- "blue"
strain.pal[names(strain.pal) %in% "JU775"] <- "green"
strain.pal[names(strain.pal) %in% "MY16"] <- "#67697C"

# Sentinel tracking
tx.sentinel.phenotypes <- toxin.sentinel.df %>%
  dplyr::filter(drug %in% gwa.metadata[which(gwa.metadata$class == "Toxicant"),]$drug) %>%
  ggplot(., mapping = aes(x = Metadata_Experiment, y = median_wormlength_um_delta_reg, fill = strain)) + 
  theme_bw() + 
  geom_boxplot(outlier.shape = NA, alpha = 0.1) + 
  stat_summary(fun = median,
               geom = "line",
               aes(group = strain, colour = strain)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = strain.pal) + 
  scale_colour_manual(values = strain.pal) + 
  facet_wrap(.~drug, scales = "free_y") + 
  labs(y = "Normalized animal length")
ggsave(tx.sentinel.phenotypes, filename = paste0("plots/sentinelphenos.", today, ".png"), width = 10, height = 7)














####################
### HERITABILITY ###
####################


H2 <- function(d){
  strain.fit <- lme4::lmer(data = d, formula = value ~ 1 + (1|strain))
  variances <- lme4::VarCorr(x = strain.fit)
  A <- as.data.frame(variances)
  Vg <- A$vcov[1]
  Ve <- A$vcov[2]
  H2 <- Vg/(Vg+Ve)
  return(H2)
}
h2 <- function(d, geno_matrix){
  pheno_strains <- unique(d$strain)
  A <- sommer::A.mat(t(geno_matrix[,colnames(geno_matrix) %in% pheno_strains]))
  
  df_y <- d %>%
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value) %>%
    dplyr::mutate(strain = as.character(strain))
  h2_res <- sommer::mmer(value~1,
                         random=~vs(strain,Gu=A), 
                         data=df_y)
  h2 <- vpredict(h2_res, h2 ~ (V1) / ( V1+V2))[[1]][1]
  h2_SE <- vpredict(h2_res, h2 ~ (V1) / ( V1+V2))[[2]][1]
  
  h2_df <- data.frame(h2) %>%
    dplyr::mutate(h2.upper = h2 + h2_SE,
                  h2.lower = h2 - h2_SE)
  return(h2_df)
}
H2.bootstrapping.calc <- function(d, nreps = 100, boot = T, geno_matrix){
  
  if(boot == T){
    # Broad-sense Heritability
    H2.point <- H2(d = d)
    h2.point <- h2(d = d, geno_matrix = geno_matrix)
    H2.boots <- list()
    for(i in 1:nreps) {
      if(i %% 10 == 0){
        print(paste0((i/nreps)*100,"%"))
      }
      #################################
      # Bootstrap within strain ##
      #################################
      nested <- d %>%
        dplyr::group_by(strain) %>%
        tidyr::nest()
      boot.strain <- list()
      for(j in 1:length(nested$strain)){
        boot.strain[[j]] <- nested$data[[j]][sample(seq(1:nrow(nested$data[[j]])),replace = T),] %>%
          dplyr::mutate(strain = nested$strain[[j]])
      }
      boot <- boot.strain %>%
        Reduce(rbind,.)
      
      ##################################
      ## Bootstrap agnostic of strain ##
      ##################################
      # boot <- d[sample(seq(1:nrow(d)), replace = T),]
      
      check <- boot %>%
        dplyr::group_by(strain) %>%
        dplyr::summarise(n())
      if(1 %in% check$`n()`){
        print("Only 1 Strain Sampled in Bootstrap - Skipping")
        next
      }
      # Broad-Sense Heritability
      H2.est <- H2(d = boot)
      H2.boots[i] <- H2.est
    }
    
    H2.boots.vec <- unlist(H2.boots)
    H2.quantiles <- c(quantile(H2.boots.vec, probs = seq(0,1,0.05)))
    H2.CI <- data.frame(H2.point, 
                        as.numeric(H2.quantiles[2]), 
                        as.numeric(H2.quantiles[21])) %>%
      `colnames<-`(c("H2.Point.Estimate","H2.5.perc","H2.95.perc"))
    
    
    
    return(list(H2.CI,H2.boots.vec,h2.point))
    
  } else {
    
    H2.point <- H2(d = d)
    # h2.point <- h2(d = d)
    H2.CI <- data.frame(H2.point, 
                        NA, 
                        NA) %>%
      `colnames<-`(c("H2.Point.Estimate","H2.5.perc","H2.95.perc"))
    return(H2.CI)
  }
  
}



all.H2s %>%
  Reduce(rbind,.) %>%
  ggplot(., mapping = aes(x = h2)) + 
  theme_bw() + 
  geom_histogram(bins = 15, fill = "#4E2A84") + 
  theme(panel.grid = element_blank()) + 
  xlim(c(0,1)) + 
  labs(x = expression(italic(h^2)),
       y = "Number of traits")

# H2.plot <- function(H2.out){
all.H2s <- list()
for(i in 1:length(all.processed.GWA.data[2:28])){
  if(is.character(all.processed.GWA.data[[i]])){
    next
  } else {
    tx <- unique(all.processed.GWA.data[[i]]$drug)
    print(tx)
    slim.phenos <- all.processed.GWA.data[[i]] %>%
      data.frame() %>%
      dplyr::select(strain, median_wormlength_um_delta_reg) %>%
      dplyr::rename(value = median_wormlength_um_delta_reg)
    H2.output <- H2.bootstrapping.calc(d = slim.phenos, geno_matrix = complete.genos)
    H2.output[[1]]$tx <- tx
    H2.output[[3]]$tx <- tx
    all.H2s[[i]] <- H2.output[[1]] %>%
      dplyr::full_join(.,H2.output[[3]]) %>%
      dplyr::rename(H2 = H2.Point.Estimate,
                    H2.lower = H2.5.perc,
                    H2.upper = H2.95.perc) %>%
      dplyr::select(contains("2"),tx)
  }
  
}
herits.df <- all.H2s %>%
  Reduce(rbind,.) %>%
  dplyr::mutate(label = if_else(h2 > 0.65 | h2 < 0.25, true = tx, false = "")) %>%
  dplyr::rename(drug = tx) %>%
  dplyr::left_join(., gwa.metadata) %>%
  dplyr::filter(class == "Toxicant") %>%
  dplyr::mutate(drug = if_else(drug == "2_4_D", true = "2,4-D", false = drug),
                drug = if_else(drug == "Nickel dichloride", true = "Nickel chloride", false = drug),
                drug = if_else(drug == "Cadmium dichloride", true = "Cadmium chloride", false = drug),
                drug = if_else(drug == "Zinc dichloride", true = "Zinc chloride", false = drug),
                
                drug = if_else(drug == "Triphenyl_phosphate 6.25", true = "Triphenyl phosphate; 6.25 μM", false = drug),
                drug = if_else(drug == "Triphenyl_phosphate 50", true = "Triphenyl phosphate; 50 μM", false = drug),
                
                drug = if_else(drug == "Paraquat_62.5", true = "Paraquat; 62.5 μM", false = drug),
                drug = if_else(drug == "Paraquat_250", true = "Paraquat; 250 μM", false = drug),
                
                drug = if_else(drug == "Silver_250", true = "Paraquat; 62.5 μM", false = drug),
                drug = if_else(drug == "Silver_250", true = "Paraquat; 250 μM", false = drug))

# H_H <- ggplot(data = herits.df, mapping = aes(x = H2, y = h2)) + 
#   theme_bw() + 
#   geom_point(aes(color = class)) + 
#   geom_errorbar(aes(ymin = h2.lower, ymax = h2.upper, color = class)) +
#   geom_errorbarh(aes(xmin = H2.lower,xmax = H2.upper, color = class)) +
#   # geom_text_repel(aes(label = label), box.padding = 1.2, max.overlaps = Inf) + 
#   geom_abline(slope = 1)  +
#   scale_colour_manual(values = c("#7E78D2","#2CA081"), name = "Compound Type") + 
#   theme(panel.grid = element_blank(),
#         legend.position = "none") + 
#   ylim(c(0,1)) + 
#   xlim(c(0,1)) + 
#   labs(x = bquote(H^2),
#        y = bquote(h^2))
# H_H
# ggsave(H_H, filename = paste("plots/HHplot", today, "png", sep = "."), width = 7, height = 7)


# H2.plot <- ggplot(data = herits.df, mapping = aes(x = reorder(drug,H2), y = H2, color = class)) + 
#   theme_bw() + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = H2.lower,ymax = H2.upper)) +
#   # geom_text_repel(aes(label = label), box.padding = 1.2, max.overlaps = Inf) + 
#   scale_colour_manual(values = c("#7E78D2","#2CA081"), name = "Compound Type") + 
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = "top") + 
#   ylim(c(0,1)) + 
#   labs(y = bquote(H^2),
#        x = "Compound")
# ggsave(H2.plot, filename = paste("plots/broadsense", today, "png", sep = "."), width = 7, height = 7)

h2.plot <- ggplot(data = herits.df, mapping = aes(x = reorder(drug,h2), y = h2)) + 
  theme_bw() + 
  geom_point() + 
  geom_errorbar(aes(ymin = h2.lower,ymax = h2.upper)) +
  # geom_text_repel(aes(label = label), box.padding = 1.2, max.overlaps = Inf) + 
  # scale_colour_manual(values = c("#7E78D2","#2CA081"), name = "Compound Type") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
        legend.position = "top") + 
  ylim(c(0,1)) + 
  labs(y = expression(italic(h^2)),
       x = "Compound")
ggsave(h2.plot, filename = paste("plots/narrowsense", today, "png", sep = "."), width = 10, height = 7)


# Diagnostics for Data Quality
# by strain
well.stats <- prelim.GWA %>%
  dplyr::select(strain, drug, Metadata_Well, Assay_Bleach, bleach, n) %>%
  # dplyr::filter(!is.na(n)) %>%
  dplyr::group_by(Assay_Bleach, strain) %>%
  dplyr::summarise(num.wells = n(),
                   perc.wells = (num.wells/228)*100, # 4 wells * 57 compounds = 228
                   mean.n = mean(n),
                   cv.n = var(n)/mean(n))

well.stats[well.stats$perc.wells <95,]


well.stat.plot <- ggplot(well.stats, mapping = aes(x = strain, y = perc.wells)) + 
  theme_bw() + 
  geom_col() + 
  facet_wrap(.~Assay_Bleach, scales = "free") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  labs(x = "Strain",
       y = "Percent of Expected Wells Retained in Final Dataset")
ggsave(well.stat.plot, filename = paste("plots/well.stat.plot.strain", today, "png", sep = "."))


# cv.n.plot <- ggplot(well.stats, mapping = aes(x = strain, y = cv.n)) + 
#   theme_bw() + 
#   geom_point() +
#   facet_wrap(.~Assay_Bleach, scales = "free_x") + 
#   ylim(c(0,20)) + 
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) + 
#   labs(x = "Strain",
#        y = "CV(# animals/well)")
# ggsave(cv.n.plot, filename = paste("plots/cv.n.plot.strain", today, "png", sep = "."), width = 9, height = 6)

well.stats %>%
  dplyr::filter(cv.n > quantile(well.stats$cv.n, probs = seq(0,1,0.05))[20]) %>%
  dplyr::group_by(Assay_Bleach) %>%
  dplyr::summarise(n())


# by compound
well.stats <- prelim.GWA %>%
  dplyr::select(strain, drug, Metadata_Well, Assay_Bleach, bleach, n) %>%
  # dplyr::filter(!is.na(n)) %>%
  dplyr::group_by(Assay_Bleach, drug) %>%
  dplyr::summarise(num.wells = n(),
                   perc.wells = (num.wells/96)*100, # 96 wells per compound per bleach
                   mean.n = mean(n),
                   cv.n = var(n)/mean(n)) %>%
  dplyr::mutate(label = if_else(drug %in% c("DMSO","Water"), true = drug, false = ""))
well.stat.plot <- ggplot(well.stats[which(!well.stats$drug %in% c("DMSO","Water")),], 
                         mapping = aes(x = drug, y = perc.wells)) + 
  theme_bw() + 
  geom_col() + 
  facet_wrap(.~Assay_Bleach, scales = "free") + 
  geom_label_repel(aes(label = label)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  labs(x = "Compound",
       y = "Percent of Expected Wells Retained in Final Dataset")
ggsave(well.stat.plot, filename = paste("plots/well.stat.plot.tx", today, "png", sep = "."))

cv.n.plot <- ggplot(well.stats[which(!well.stats$drug %in% c("DMSO","Water")),], mapping = aes(x = drug, y = cv.n)) + 
  theme_bw() + 
  geom_point() +
  facet_wrap(.~Assay_Bleach, scales = "free_x") + 
  # ylim(c(0,20)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  labs(x = "Compound",
       y = "CV(# animals/well)")
ggsave(cv.n.plot, filename = paste("plots/cv.n.plot.tx", today, "png", sep = "."), width = 9, height = 6)




