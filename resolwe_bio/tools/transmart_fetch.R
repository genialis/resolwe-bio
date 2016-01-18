#!/usr/bin/Rscript
require("transmartRClient")
require('argparse')

parser = ArgumentParser(description='Fetch data from tranSMART instance.')
parser$add_argument('--URL', default = "http://137.117.169.225:8080/transmart", help='tranSMART instance URL')
parser$add_argument('--annConceptLinks', default = "", help='Annotations concept links')
parser$add_argument('--expsConceptLinks', default = "", help='Expressions concept links')
parser$add_argument('--token', default = "", help='Auth token')
parser$add_argument('--projection', default='log_intensity', help='Name of the data projection to fetch (log_intensity, all_data, default_real_projection, ...')
parser$add_argument('--outA', help='Output annotation file name')
parser$add_argument('--outE', help='Output expression file name')
parser$add_argument('--outT', help='Output tree file name')

args = parser$parse_args(commandArgs(trailingOnly=TRUE))

# TODO: Why is total NA and why do current and total have 2 values???
.downloadcallback <- function() {
    start <- function(.total) cat("Retrieving data...\n")
    update <- function(current, total) {
        if (current[2] > 0 && !is.na(total[2]) && total[2] > 0) {
            # limit the progress to 0.95 (5 % left for parsing)
            cat('{"proc.progress":', current[2] / total[2] * 0.95, '}\n')
        }
    }
    end <- function() cat("Download complete.\n")
    environment()
}

connectToTransmart(args$URL, .access.token = args$token)

# get annotations
if (args$annConceptLinks != '') {
    ann <- readChar(args$annConceptLinks, file.info(args$annConceptLinks)$size)
    ann <- gsub("[\r\n]", "", ann)
    links_new <- c(unlist(strsplit(ann, ';')))

    observations <- getObservations(concept.links = links_new, as.data.frame = T)

    final <- data.frame(cbind(observations$subjectInfo$subject.inTrialId, observations$observations))
    final <- final[, !names(final) == 'subject.id']
    colnames(final)[1] <- 'ID'
    final <- as.data.frame(lapply(final, as.character), stringsAsFactors = FALSE)

    empty <- rep("", ncol(final))
    final <- rbind(empty, final)
    final <- rbind(empty, final)

    write.table(final, file = args$outA, quote = FALSE, sep = "\t", row.names = F, na = "")

    # get study tree
    con <- c()
    stu <- list()
    for (study in c(unlist(strsplit(ann, ';')))) {
        study <- c(strsplit(as.character(study), '/'))
        studyId <- study[[1]][3]

        if (!(studyId %in% stu)) {
            studyConcepts = getConcepts(studyId)
            con <- c(con, studyConcepts$fullName)
            stu <- c(stu, studyId)
        }
    }

    write.table(con, file = args$outT, quote = FALSE, sep = "\t", row.names = F, col.names = F)
}

# get expressions
if (args$expsConceptLinks != '') {
    dataDownloaded <- getHighdimData(concept.link = args$expsConceptLinks, projection=args$projection, progress.download = .downloadcallback())

    data = dataDownloaded[[1]]
    expression_data = data[,-c(1:7)]
    rownames(expression_data) = make.names(data$patientId, unique=TRUE)
    #rownames(expression_data) = data$patientId
    write.table(t(expression_data), file = args$outE, quote = FALSE, sep = "\t", col.names = NA)
}
