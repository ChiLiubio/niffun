#' @title
#' Create an R6 object for functional prediction of nitrogen fixers based on nifH sequences.
#'
#' @description
#' This class is a wrapper for a series of functional analysis on diazotrophic communities based on the nifH sequences.
#'
#' @export
trans_niffun <- R6::R6Class(classname = "trans_niffun",
	public = list(
		#' @description
		#' Create the trans_niffun object.
		#' 
		#' @param dataset default NULL; the object of microtable Class.
		#' @examples
		#' \donttest{
		#' data(dataset_nifH)
		#' t1 <- trans_niffun$new(dataset = dataset)
		#' }
		initialize = function(dataset = NULL){
			if(is.null(dataset)){
				stop("Please provide dataset !")
			}
			if(is.null(dataset$rep_fasta)){
				stop("The rep_fasta is missing in your dataset object! The nifH fasta file is necessary! Please use help(microtable) to see the rep_fasta description!")
			}
			self$dataset <- clone(dataset)
		},
		#' @description
		#' Taxonomic assignment with blast.
		#' 
		#' @param blast_tool_path default NULL; the folder path, e.g. ncbi-blast-2.5.0+/bin ; blast tools folder downloaded from 
		#'   "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+"  ; e.g. ncbi-blast-2.5.0+-x64-win64.tar.gz  for windows system; 
		#'   if blast_tool_path is NULL, search the tools in the environmental path variable.
		#' @param path_to_temp_folder default NULL; The temporary folder to store the logfile, intermediate file and result files; if NULL, 
		#' 	 use the default temporary in the computer.
		#' @param evalue default 1e-5; the E value threshold used in blast.
		#' @param max_target_seqs default 50; Maximum number of aligned sequences to keep.
		#' @param num_threads default 1; Number of threads (CPUs) to use in the BLAST search.
		#' @param other_blast_parameters default ""; other parameters provided to BLAST, such as "-perc_identity 80".
		#' @return res_tax_table stored in the object. 
		#' @examples
		#' \donttest{
		#' t1$cal_blast()
		#' }
		cal_blast = function(
			blast_tool_path = NULL,
			path_to_temp_folder = NULL,
			evalue = 1e-5,
			max_target_seqs = 50,
			num_threads = 1,
			other_blast_parameters = ""
			){
			if(is.null(path_to_temp_folder)){
				path_to_temp_folder <- tempdir()
				message("The intermediate files are saved in the temporary directory --- ", path_to_temp_folder)
			}
			if(!dir.exists(path_to_temp_folder)){
				stop(paste0("Temporay folder--", path_to_temp_folder, " is not existed! Please check it!"))
			}
			
			# judge input type: DNA or protein
			unique_char <- self$dataset$rep_fasta[[1]] %>% 
				as.character %>% 
				strsplit(split = "") %>% 
				unlist %>%
				unique
			if(all(unique_char %in% c("G", "A", "C", "T", "N"))){
				seq_type <- "DNA"
			}else{
				seq_type <- "protein"
			}
		
			# check whether blast tool is available
			blast_soft <- ifelse(seq_type == "DNA", "blastn", "blastp")
			if(is.null(blast_tool_path)){
				blast_bin <- blast_soft
			}else{
				if (tolower(Sys.info()[["sysname"]]) != "windows"){
					blast_bin <- file.path(blast_tool_path, blast_soft)
					system(paste("chmod +x", blast_bin))
				}else{
					blast_bin <- file.path(blast_tool_path, paste0(blast_soft, ".exe"))
				}
			}
			res <- system(command = paste(blast_bin, "-help"), intern = T)
			if(length(res) == 0){
				blast_bin <- blast_soft
				res <- system(command = paste(blast_bin, "-help"), intern = T)
				if(length(res) == 0){
					stop(paste0(blast_bin, " not found! Please check the file path!"))
				}
			}
			# copy database to the temp dir
			ref_db_file <- ifelse(seq_type == "DNA", "nifH_DNA.fasta", "nifH_protein.fasta")
			path_from_ref_db <- system.file("extdata", ref_db_file, package = "niffun")
			path_to_ref_db <- file.path(path_to_temp_folder, ref_db_file)
			file.copy(from = path_from_ref_db, to = path_to_ref_db)
			# makeblastdb
			if(is.null(blast_tool_path)){
				makeblastdb_bin <- "makeblastdb"
			}else{
				makeblastdb_bin <- file.path(blast_tool_path, "makeblastdb")
				if (tolower(Sys.info()[["sysname"]]) != "windows"){
					system(paste("chmod +x", makeblastdb_bin))
				}else{
					makeblastdb_bin <- file.path(blast_tool_path, "makeblastdb.exe")
				}
			}
			res <- system(command = paste(makeblastdb_bin, "-help"), intern = T)
			if(length(res) == 0){
				makeblastdb_bin <- "makeblastdb"
				res <- system(command = paste(makeblastdb_bin, "-help"), intern = T)
				if(length(res) == 0){
					stop(paste0(makeblastdb_bin, " not found! Please check the file path!"))
				}
			}
			# use refernence db
			message('Generate blast reference database ...')
			makedb_type <- ifelse(seq_type == "DNA", "nucl", "prot")
			cmd <- paste(makeblastdb_bin, '-dbtype', makedb_type, '-in', path_to_ref_db)
			if(tolower(Sys.info()[["sysname"]]) == "windows"){
				system(cmd, show.output.on.console = F)
			}else{
				system(cmd, ignore.stdout = T, ignore.stderr = T)
			}
			message('Reference database finished ...')
			# write the fasta file
			rep_fasta_path <- file.path(path_to_temp_folder, "rep_fasta.tmp.fasta")
			if(inherits(self$dataset$rep_fasta, "list")){
				seqinr::write.fasta(self$dataset$rep_fasta, names = names(self$dataset$rep_fasta), file.out = rep_fasta_path)
			}else{
				if(inherits(self$dataset$rep_fasta, "DNAStringSet")){
					Biostrings::writeXStringSet(x = self$dataset$rep_fasta, filepath = rep_fasta_path)
				}else{
					stop("Unknown fasta format! Must be either list (from read.fasta of seqinr package) or DNAStringSet (from readDNAStringSet of Biostrings package) !")
				}
			}
			
			res_blast_path <- file.path(path_to_temp_folder, 'ref_blast.txt')
			message('Blast start ...')
			cmd <- paste(blast_bin, '-db', path_to_ref_db, '-query', rep_fasta_path, '-evalue', evalue, '-max_target_seqs', max_target_seqs, 
				'-outfmt 6 -out', res_blast_path, '-num_threads', num_threads, other_blast_parameters)
			if (tolower(Sys.info()[["sysname"]]) == "windows"){
				system(cmd, show.output.on.console = F)
			}else{
				system(cmd, ignore.stdout = T, ignore.stderr = T)
			}
			message('Blast finished ...')
			# read the blast file
			ref_blast_result <- read.delim(res_blast_path, header = F)
			# taxonomic assignment
			data("nifH_seqtag_genometag", envir = environment())
			use_join_name <- ifelse(seq_type == "DNA", "DNA_name", "protein_name")
			total_info <- left_join(ref_blast_result, nifH_seqtag_genometag, by = c("V2" = use_join_name))
			total_info %<>% .[!duplicated(.[, 1]), ]

			total_info$Species[total_info[, 3] >= 88.1 & total_info[, 3] < 91.9] <- ""
			select_rows <- total_info[, 3] >= 75 & total_info[, 3] < 88.1
			if(any(select_rows)){
				total_info$Genus[select_rows] <- ""
				total_info$Species[select_rows] <- ""
			}
			select_rows <- total_info[, 3] < 75
			if(any(select_rows)){
				total_info$Family[select_rows] <- ""
				total_info$Genus[select_rows] <- ""
				total_info$Species[select_rows] <- ""
			}

			tax_assignment <- total_info[, c("V1", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
			rownames(tax_assignment) <- tax_assignment[, 1]
			tax_assignment %<>% .[, -1]
			self$res_blast_tax <- total_info
			message('Blast raw result is stored in object$res_tax_table ...')
			self$res_tax_table <- tax_assignment
			message('Taxonomic table is stored in object$res_tax_table ...')
			self$path_to_temp_folder <- path_to_temp_folder
		},
		#' @description
		#' Predict KEGG orthology and pathway abundance.
		#' 
		#' @param min_identity_to_reference default 95; the sequences identity threshold used for finding the nearest species genome.
		#' @return res_KEGG_KO and res_KEGG_pathway stored in the object. 
		#' @examples
		#' \donttest{
		#' t1$cal_pathway()
		#' }
		cal_pathway = function(min_identity_to_reference = 95){
			if(is.null(self$path_to_temp_folder)){
				path_to_temp_folder <- tempdir()
				message("The intermediate files are saved in the temporary directory --- ", path_to_temp_folder)
			}else{
				path_to_temp_folder <- self$path_to_temp_folder
			}
			path_to_log_file <- file.path(path_to_temp_folder, 'logfile.txt')
			if(min_identity_to_reference < 88){
				warning("Minimum identity of less than 88% (minimum similarity threshold of Genera) will likly results in inaccurate predictions!")
			}
			message(paste0("Using minimum identity cutoff of ", min_identity_to_reference, "% to nearest neighbor"))
			ref_blast_result <- self$res_blast_tax
			ref_blast_result_reduced <- ref_blast_result[which(ref_blast_result$V3 >= min_identity_to_reference), c("V1", "V2", "genome_tag")]

			# Reading and filtering the otu table
			otu_table <- self$dataset$otu_table %>% tibble::rownames_to_column()
			# for the taxa mapping
			raw_otu_table_reduced <- merge(x = ref_blast_result_reduced, y = otu_table, by.x = "V1", by.y = "rowname")
			otu_table_reduced <- raw_otu_table_reduced[, -c(1:2)]
			# for the calculation
			otu_table_reduced_aggregated <- aggregate(x = otu_table_reduced[, -1, drop = FALSE], by = list(otu_table_reduced[,1]), sum)

			# Write unknown fraction to log file
			unknown_fraction1 <- round(1 - colSums(ifelse(otu_table_reduced[, -1, drop = FALSE] > 0, 1, 0)) / colSums(ifelse(otu_table[, -1, drop = FALSE] > 0, 1, 0)), digits = 5) %>%
				as.data.frame
			write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
			write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
			unknown_fraction2 = round(1 - colSums(otu_table_reduced_aggregated[, -1, drop = FALSE]) / colSums(otu_table[, -1, drop = FALSE]), digits = 5) %>%
				as.data.frame
			write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
			write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)

			# Generate reference profile
			data("mapped_KO", envir = environment())
			reference_profile <- mapped_KO %>% t %>% as.data.frame %>% .[otu_table_reduced_aggregated$Group.1, ]

			# all the required KEGG files are stored in Tax4Fun2_KEGG Rdata
			data("Tax4Fun2_KEGG", package = "microeco", envir=environment())
			ko_list <- Tax4Fun2_KEGG$ko_list

			map_reference_profile <- reference_profile %>% tibble::rownames_to_column(var = "id")
			res_reference_profile <- dplyr::left_join(raw_otu_table_reduced[, c(1:3)], map_reference_profile, by = c("genome_tag" = "id"))
			colnames(res_reference_profile)[1:2] <- c("Taxa", "seq_tag")
			write.table(res_reference_profile, file.path(path_to_temp_folder, 'res_reference_profile.tsv'), sep = "\t")
			message("Reference profile file is saved in ", file.path(path_to_temp_folder, "res_reference_profile.tsv"), " ...")

			# Calculate functional profiles sample-wise
			message('Generating functional profile for samples ...')
			functional_prediction <- NULL
			for(sample_num in 2:ncol(otu_table_reduced_aggregated)){
				functional_prediction_sample <- reference_profile * as.numeric(otu_table_reduced_aggregated[, sample_num])
				functional_prediction_sample <- colMeans(functional_prediction_sample)
				functional_prediction_sample <- functional_prediction_sample / sum(functional_prediction_sample)
				if(is.na(sum(functional_prediction_sample))){
					functional_prediction_sample[1:length(functional_prediction_sample)] <- 0
				}
				functional_prediction <- cbind(functional_prediction, functional_prediction_sample)
			}

			colnames(functional_prediction) <- names(otu_table)[-1]
			functional_prediction_final <- data.frame(KO = rownames(functional_prediction), functional_prediction, 
				description = ko_list$description[match(rownames(functional_prediction), ko_list$ko)])
			if(ncol(functional_prediction) >= 2) keep <- which(rowSums(functional_prediction) > 0)
			if(ncol(functional_prediction) == 1) keep <- which(functional_prediction > 0)
			if (length(keep) == 0){
				stop("No functional prediction possible!\nEither no nearest neighbor found or your table is empty!")
			}
			functional_prediction_final <- functional_prediction_final[keep, ]
			write.table(x = functional_prediction_final, file = file.path(path_to_temp_folder, 'KEGG_KO.txt'), append = F, quote = F, 
				sep = "\t", row.names = F, col.names = T)
			self$res_KEGG_KO <- functional_prediction_final
			message('Result KO abundance is stored in object$res_KEGG_KO ...')

			# Converting the KO profile to a profile of KEGG pathways
			message('Converting functions to pathways')
			ko2ptw <- Tax4Fun2_KEGG$ko2ptw
			ko2ptw$nrow <- Tax4Fun2_KEGG$ko_list[ko2ptw$nrow, "ko"]
			ko2ptw %<>% .[.$nrow %in% rownames(functional_prediction), ]
			pathway_prediction <- aggregate(x = functional_prediction[ko2ptw$nrow, ], by = list(ko2ptw$ptw), sum)

			col_sums <- colSums(pathway_prediction[, -1, drop = FALSE])
			col_sums[col_sums == 0] <- 1
			pathway_prediction[, -1] <- t(t(pathway_prediction[, -1, drop = FALSE]) / col_sums)
			keep <- which(rowSums(pathway_prediction[, -1, drop = FALSE]) > 0)

			if(sum(pathway_prediction[, -1]) == 0) stop("Conversion to pathway failed! No pathway abundance obtained!")
			pathway_prediction %<>% .[keep, , drop = FALSE]
			rownames(pathway_prediction) <- pathway_prediction[, 1]
			pathway_prediction <- pathway_prediction[, -1, drop = FALSE]

			self$res_KEGG_pathway <- pathway_prediction
			message('KEGG pathway abundance table is stored in object$res_KEGG_pathway ...')			
			ptw_desc <- Tax4Fun2_KEGG$ptw_desc
			pathway_prediction_final <- data.frame(pathway_prediction, ptw_desc[rownames(pathway_prediction), ])
			pathway_prediction_final <- data.frame(pathway = rownames(pathway_prediction_final), pathway_prediction_final)
			write.table(x = pathway_prediction_final, file = file.path(path_to_temp_folder, 'KEGG_pathway_prediction.txt'), 
				append = F, quote = F, sep = "\t", row.names = F, col.names = T)

		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)






