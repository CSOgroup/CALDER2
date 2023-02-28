# Testing using Calder with a .cool file as input
sanitize_chroms <- function(chroms){
    res <- lapply(chroms, function(x){
        if(startsWith(x, "chr")){
            return(substring(x, 4))
        } else{
            return(x)
        }
    })
    return(res)
}


handle_input_cool <- function(input, 
                              outpath, 
                              bin_size=50000, 
                              genome="hg38",  
                              nproc=10,
                              chroms_to_remove = c("MT", "M", 'chrMT', 'chrM', 'Y', 'chrY')){

    intermediate_data_dir = file.path(outpath, "intermediate_data")
    dir.create(intermediate_data_dir, recursive=TRUE, showWarnings=FALSE)

    system(paste0("cooler dump --table chroms --out ", 
                  file.path(intermediate_data_dir, "chroms.txt"), 
                  " --header ", 
                  input))
    chroms <- read.table(file.path(intermediate_data_dir, "chroms.txt"), sep="\t", header=TRUE)
    chroms <- chroms[!(chroms$name %in% chroms_to_remove), "name"]

    dump_paths <- list()
    for(chrom in chroms){
        cat(paste0("[Pre-processing] Dumping ", chrom, "\n"))
        chrom_dump_path <- file.path(intermediate_data_dir, paste0(chrom, "_dump.txt"))
        dump_paths <- c(dump_paths, chrom_dump_path)
        if(! file.exists(chrom_dump_path)){
            system(paste0("cooler dump --table pixels --range ", 
                          chrom, 
                          " --join --balanced ",
                          input,
                          " | cut -f2,5,8 | awk '{if ($3) print;}' > ",
                          chrom_dump_path))
        }
    }

    chroms <- sanitize_chroms(chroms)
    names(dump_paths) <- chroms

    CALDER(contact_file_dump=dump_paths, 
           chrs=chroms, 
           bin_size=bin_size,
           genome=genome,
           save_dir=outpath,,
           single_binsize_only=TRUE,
           save_intermediate_data=TRUE,
           n_cores=nproc,
           sub_domains=TRUE)
    file.remove(file.path(testthat::test_path(), "total_execution.time"))
}

test_that("CALDER works with cool files", {
    input_cool_path <- file.path(testthat::test_path("data"), "test.cool")
    output_path <- testthat::test_path("test-main-cool-out")
    handle_input_cool(input_cool_path, output_path)
    expect_snapshot_file(file.path(output_path, "sub_compartments", "all_sub_compartments.bed"), name = "TestCool_all_sub_compartments.bed")
    expect_snapshot_file(file.path(output_path, "sub_compartments", "all_sub_compartments.tsv"), name = "TestCool_all_sub_compartments.tsv")
    expect_snapshot_file(file.path(output_path, "sub_domains", "all_nested_boundaries.bed"), name = "TestCool_all_nested_boundaries.bed")
    unlink(output_path, recursive=TRUE)
})


test_that("CALDER works with dumps", {

    chrs = c(21:22)

    ## demo contact matrices in dump format
    contact_file_dump = as.list(system.file("extdata", sprintf("mat_chr%s_10kb_ob.txt.gz", chrs),
                package='CALDER'))
    names(contact_file_dump) = chrs
    output_path <- testthat::test_path("test-main-dump-out")
    ## Run CALDER to compute compartments but not nested sub-domains
    CALDER(contact_file_dump=contact_file_dump, 
            chrs=chrs, 
            bin_size=10E3,
            genome='hg19',
            save_dir=output_path,
            save_intermediate_data=FALSE,
            n_cores=2,
            sub_domains=FALSE)
    expect_snapshot_file(file.path(output_path, "sub_compartments", "all_sub_compartments.bed"), name = "TestDump_all_sub_compartments.bed")
    expect_snapshot_file(file.path(output_path, "sub_compartments", "all_sub_compartments.tsv"), name = "TestDump_all_sub_compartments.tsv")
    unlink(output_path, recursive=TRUE)
})

