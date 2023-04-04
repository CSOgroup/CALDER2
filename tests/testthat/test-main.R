
CALDER_CLI=file.path(system.file("scripts", package="CALDER"), "calder")


test_that("CALDER works with cool files", {
    input_cool_path <- file.path(testthat::test_path("data"), "test.cool")
    output_path <- testthat::test_path("test-main-cool-out")

    CMD = paste0(CALDER_CLI, " -i ", input_cool_path, " -t cool -g hg38 -o ", output_path)
    system(CMD)

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

test_that("CALDER works with custom feature track", {
    input_cool_path <- file.path(testthat::test_path("data"), "test.cool")
    feature_track_path <- file.path(testthat::test_path("data"), "test_gene_coverage.bed")
    output_path <- testthat::test_path("test-main-featuretrack-out")

    CMD = paste0(CALDER_CLI, " -i ", input_cool_path, " -t cool -g hg38 -o ", output_path, " -f ", feature_track_path)
    system(CMD)

    expect_snapshot_file(file.path(output_path, "sub_compartments", "all_sub_compartments.bed"), name = "TestFeatureTrack_all_sub_compartments.bed")
    expect_snapshot_file(file.path(output_path, "sub_compartments", "all_sub_compartments.tsv"), name = "TestFeatureTrack_all_sub_compartments.tsv")
    expect_snapshot_file(file.path(output_path, "sub_domains", "all_nested_boundaries.bed"), name = "TestFeatureTrack_all_nested_boundaries.bed")
    unlink(output_path, recursive=TRUE)
})