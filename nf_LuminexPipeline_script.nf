nextflow.enable.dsl=2

params.aref = "/home/jesse-build/nf/aref/lum_analyte_ref.RData"
params.data_dir = '~/Documents/datasets'
params.tech_reps = 1
params.instrument_names = 'bp, mp'
//params.instrument_name2 = "mp"

process CONFIG{
  
  input:
  val data_dir
  val aref
  val tech_reps
  
  script:
  """
  #!/usr/bin/env Rscript
  
  setwd("${projectDir}") 
  
  arefs <- "$aref"

  LuminexPipeline::pipeline_config(wd = "${projectDir}", dd = "$data_dir", aref = arefs, tech_reps = $tech_reps) 
  
  """
}
  process DATA_IMPORT{
  
  input:
  val data_dir
  
  output:
  val "${projectDir}/rds/dta_import.rds"
  
  script:
  """
  #!/usr/bin/env Rscript

  setwd("${projectDir}")
  
  LuminexPipeline::data_import("$data_dir") 
  """
  }

  process FILENAME_SEPARATE{
    input:
    val "data_out"
    val instrument_names
    //val instrument_name2

    output:
    val "${projectDir}/rds/dta_separate.rds"

    script:
    """
    #!/usr/bin/env Rscript

    setwd("${projectDir}")
    
    dta <- readRDS("${data_out}")

    instrument_names <- c("$instrument_names")

    LuminexPipeline::filename_separate(data = dta, instrument_names = instrument_names)
    """
  }
  
  process DATA_CLEAN{
  
  input:
  val "data_out"

  output:
  val "${projectDir}/rds/dta_colnames_clean.rds"

  script:
  """
  #!/usr/bin/env Rscript

  setwd("${projectDir}")  

  dta1 <- readRDS("${data_out}")

  LuminexPipeline::colnames_clean(dta1)

  dta2 <- readRDS("${data_out}")
  """
}

process ANALYTE_FIX{
  
  input:
  val "data_out"
  val aref

  output:
  val "${projectDir}/rds/dta_analyte_ref.rds"
  
  script:
  """
  #!/usr/bin/env Rscript

  setwd("${projectDir}")

  arefs <- "$aref"

  attach("$aref")

  dta2 <- readRDS("${data_out}")
  
  LuminexPipeline::analyte_names_fix2(arefs, dta2)
  """
}

process DATA_SAVE {
  
  input:
  val "data_out"

  output:
  val "${projectDir}/rds/dta_raw.rds"
  val "${projectDir}/rds/dta_symbol_remove.rds"

  script:
  """
  #!/usr/bin/env Rscript

  setwd("${projectDir}")

  dta3 <- readRDS("${data_out}")

  LuminexPipeline::save_raw(dta3)

  LuminexPipeline::symbols_remove(dta3) 
  """
}

process DATA_SPLIT {
  input:
  val "data_out"

  output:
  val "${projectDir}/rds/dta_list.rds"

  script:
  """
  #!/usr/bin/env Rscript

  setwd("${projectDir}")

  dta4 <- readRDS("${data_out}")

  LuminexPipeline::data_split(dta = dta4)
  """

}

//data_dir = Channel.value('~/Documents/datasets')


workflow{

  def analyte_ref_ch = Channel.value(params.aref)
  def data_dir_ch = Channel.value(params.data_dir)
  def data_out = Channel.value("${projectDir}/rds/dta.rds")
  def tech_rep_ch = Channel.value(params.tech_reps)
  def ins_names = Channel.value(params.instrument_names)
  //def ins_name2 = Channel.value(params.instrument_name2)

  CONFIG(data_dir_ch, analyte_ref_ch, tech_rep_ch)
  DATA_IMPORT(data_dir_ch)
  FILENAME_SEPARATE(DATA_IMPORT.out, ins_names)
  DATA_CLEAN(FILENAME_SEPARATE.out)
  ANALYTE_FIX(DATA_CLEAN.out, analyte_ref_ch)
  DATA_SAVE(ANALYTE_FIX.out)
  DATA_SPLIT(DATA_SAVE.out[1])
}
