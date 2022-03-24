#' @title
#' Create an R6 object for phenotypic and functional prediction of nitrogen fixers.
#'
#' @description
#' This class is a wrapper for a series of phenotypic and functional analysis on diazotrophic communities based on the taxonomy.
#'
#' @export
trans_nifFun <- R6::R6Class(classname = "trans_nifFun",
	public = list(
		#' @description
		#' Create the trans_nifFun object.
		#' 
		#' @param dataset default NULL; the object of microtable Class.
		#' @examples
		#' \donttest{
		#' data(dataset_nifH)
		#' t1 <- trans_nifFun$new(dataset = dataset)
		#' }
		initialize = function(dataset = NULL){
			if(is.null(dataset)){
				stop("Please provide dataset !")
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
			
			if(is.null(self$dataset$rep_fasta)){
				stop("The rep_fasta is missing in your dataset object! The fasta file is necessary in Tax4Fun2! Use help(microtable) to see the rep_fasta description!")
			}else{
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
			path_from_ref_db <- system.file("extdata", ref_db_file, package="nifFun")
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
			if (tolower(Sys.info()[["sysname"]]) == "windows"){
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
			data("nifH_seqtag_genometag_taxonomy", envir = environment())
			use_join_name <- ifelse(seq_type == "DNA", "DNA_name", "protein_name")
			total_info <- left_join(ref_blast_result, nifH_seqtag_genometag_taxonomy, by = c("V2" = use_join_name))
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
			self$res_tax_table <- tax_assignment
			message('Taxonomic table is stored in object$tax_table ...')
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)












