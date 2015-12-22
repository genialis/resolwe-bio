#!/usr/bin/Rscript

require("transmartRClient")
require('argparse')

parser = ArgumentParser(description='Fetch data from tranSMART instance.')
parser$add_argument('--URL', default = "http://52.48.28.218:8080/transmart", help='tranSMART instance URL')
parser$add_argument('--annConceptLinks', default = "", help='Annotations concept links')
parser$add_argument('--expsConceptLinks', default = "", help='Expressions concept links')
parser$add_argument('--token', default = "", help='Auth token')
parser$add_argument('--projection', default='default_real_projection', help='Name of the data projection to fetch (log_intensity, all_data, default_real_projection, ...')
parser$add_argument('--out', help='Output file name')

args = parser$parse_args(commandArgs(trailingOnly=TRUE))

connectToTransmart(args$URL, .access.token = args$token)

# get annotations -----------------------------------------------------
if (args$annConceptLinks != '') {
  links_new <- c(unlist(strsplit(args$annConceptLinks, ';')))

  observations <- getObservations(concept.links = links_new, as.data.frame = T)
  
  final <- data.frame(cbind(observations$subjectInfo$subject.inTrialId, observations$observations))
  final <- final[, !names(final) == 'subject.id']
  colnames(final)[1] <- 'ID'
  final <- as.data.frame(lapply(final, as.character), stringsAsFactors = FALSE)
  
  empty <- rep("", ncol(final))
  final <- rbind(empty, final)
  final <- rbind(empty, final)
  
  write.table(final, file = args$out, quote = FALSE, sep = "\t", row.names = F, na = "")
}

# get expressions ---------------------------------------------------

if (args$expsConceptLinks != '') {
  dataDownloaded <- getHighdimData(concept.link = args$expsConceptLinks, projection=args$projection)
  data = dataDownloaded[[1]]
  expression_data = data[,-c(1:7)]
  rownames(expression_data) = make.names(data$patientId, unique=TRUE)
  #rownames(expression_data) = data$patientId
  write.table(t(expression_data), file = args$out, quote = FALSE, sep = "\t", col.names = NA)
}
