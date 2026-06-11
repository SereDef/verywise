fix_stacks <- function(res_dir) {
  all_stacks <-list.files(res_dir, pattern = "^stack_names.txt$", recursive = TRUE, full.names = TRUE)
  
  for (stack_file in all_stacks) {
    tab <- utils::read.table(stack_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    tab$stack_name <- sort(tab$stack_name)

    utils::write.table(tab, stack_file, sep = "\t", row.names = FALSE)
  }
  
}

tab <- fix_stacks('/Users/Serena/Desktop/Packages/verywise/tests/testthat/fixtures/')

sort(tab$stack_name)
