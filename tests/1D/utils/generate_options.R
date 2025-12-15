generate_options <- function(name_main_test, path_queue) {
  ## create the directory if it does not exist
  mkdir(c(path_queue))
  
  ## name test suite
  test_suite <- "1D"
  
  ## generate the options json files
  switch(name_main_test,
         test0 = {
           
           # Test 0:
           
           ## options grid
           n_stat_units_vect <- c(100)
           NSR_X_lc_vect <- 0
           NSR_Y_vect <- 0
           
           ## fixed options
           name_mesh <- "unit_segment"
           name_mesh_short <- "us"
           n_nodes <- 51
           n_locs <- 51
           locs_eq_nodes <- TRUE
           n_comp <- 4
           
           n_reps <- 5
           
           ## assembly jsons           
           for (n_stat_units in n_stat_units_vect) {
             for (NSR_X_lc_vect in NSR_X_lc_vect) {
               for (NSR_Y in NSR_Y_vect) {
                 name_test <- paste(
                   test_suite, name_main_test, name_mesh_short,
                   sprintf("%0*d", 4, n_nodes),
                   sprintf("%0*d", 4, n_locs),
                   sprintf("%0*d", 4, n_stat_units),
                   NSR_X_lc_vect, NSR_Y, n_comp,
                   sep = "_"
                 )
                 json_list <- list(
                   test = list(
                     name_test = name_test
                   ),
                   mesh = list(
                     name_mesh = name_mesh
                   ),
                   dimensions = list(
                     n_nodes = n_nodes,
                     n_locs = n_locs,
                     n_stat_units = n_stat_units,
                     n_comp = n_comp,
                     n_reps = n_reps
                   ),
                   data = list(
                     locs_eq_nodes = locs_eq_nodes
                   ),
                   noise = list(
                     NSR_X_lc_vect = NSR_X_lc_vect,
                     NSR_Y = NSR_Y
                   ),
                   COMPETITORS = FALSE
                 )
                 write_json(json_list, path = paste(path_queue, name_test, ".json", sep = ""))
               }
             }
           }
         },
         test1 = {
           
           # Test 1:
           
           ## options grid
           n_stat_units_vect <- c(50, 100, 200)
           NSR_X_lc_vect <- seq(0, 3, length = 4)
           NSR_Y_vect <- c(0, 0.01, 0.1)
           
           ## fixed options
           name_mesh <- "unit_segment"
           name_mesh_short <- "us"
           n_nodes <- 51
           n_locs <- 200
           locs_eq_nodes <- FALSE
           n_comp <- 4
           
           n_reps <- 30
           
           ## assembly jsons           
           for (n_stat_units in n_stat_units_vect) {
             for (NSR_X_lc in NSR_X_lc_vect) {
               for (NSR_Y in NSR_Y_vect) {
                 name_test <- paste(
                   test_suite, name_main_test, name_mesh_short,
                   sprintf("%0*d", 4, n_nodes),
                   sprintf("%0*d", 4, n_locs),
                   sprintf("%0*d", 4, n_stat_units),
                   NSR_X_lc, NSR_Y, n_comp,
                   sep = "_"
                 )
                 json_list <- list(
                   test = list(
                     name_test = name_test
                   ),
                   mesh = list(
                     name_mesh = name_mesh
                   ),
                   dimensions = list(
                     n_nodes = n_nodes,
                     n_locs = n_locs,
                     n_stat_units = n_stat_units,
                     n_comp = n_comp,
                     n_reps = n_reps
                   ),
                   data = list(
                     locs_eq_nodes = locs_eq_nodes
                   ),
                   noise = list(
                     NSR_X_lc_vect = NSR_X_lc,
                     NSR_Y = NSR_Y
                   ),
                   COMPETITORS = FALSE
                 )
                 write_json(json_list, path = paste(path_queue, name_test, ".json", sep = ""))
               }
             }
           }
         },
         test2 = {
           
           # Test 2:
           
           ## options grid
           n_locs_vect <- c(50, 100, 200)
           NSR_X_lc_vect <- seq(0, 3, length = 4)
           NSR_Y_vect <- c(0, 0.01, 0.1)
           
           ## fixed options
           name_mesh <- "unit_segment"
           name_mesh_short <- "us"
           n_nodes <- 51
           n_stat_units <- 100
           locs_eq_nodes <- FALSE
           n_comp <- 4
           
           n_reps <- 30
           
           ## assembly jsons           
           for (n_locs in n_locs_vect) {
             for (NSR_X_lc in NSR_X_lc_vect) {
               for (NSR_Y in NSR_Y_vect) {
                 name_test <- paste(
                   test_suite, name_main_test, name_mesh_short,
                   sprintf("%0*d", 4, n_nodes),
                   sprintf("%0*d", 4, n_locs),
                   sprintf("%0*d", 4, n_stat_units),
                   NSR_X_lc, NSR_Y, n_comp,
                   sep = "_"
                 )
                 json_list <- list(
                   test = list(
                     name_test = name_test
                   ),
                   mesh = list(
                     name_mesh = name_mesh
                   ),
                   dimensions = list(
                     n_nodes = n_nodes,
                     n_locs = n_locs,
                     n_stat_units = n_stat_units,
                     n_comp = n_comp,
                     n_reps = n_reps
                   ),
                   data = list(
                     locs_eq_nodes = locs_eq_nodes
                   ),
                   noise = list(
                     NSR_X_lc_vect = NSR_X_lc,
                     NSR_Y = NSR_Y
                   ),
                   COMPETITORS = FALSE
                 )
                 write_json(json_list, path = paste(path_queue, name_test, ".json", sep = ""))
               }
             }
           }
         },
         test3 = {
           
           # Test 3:
           
           ## options grid
           n_nodes_vect <- c(5, 11, 21, 51, 101)
           NSR_X_lc_vect <- seq(0, 3, length = 4)
           NSR_Y_vect <- c(0, 0.01, 0.1)
           
           ## fixed options
           name_mesh <- "unit_segment"
           name_mesh_short <- "us"
           n_locs <- 200
           n_stat_units <- 100
           locs_eq_nodes <- FALSE
           n_comp <- 4
           
           n_reps <- 30
           
           ## assembly jsons           
           for (n_nodes in n_nodes_vect) {
             for (NSR_X_lc in NSR_X_lc_vect) {
               for (NSR_Y in NSR_Y_vect) {
                 name_test <- paste(
                   test_suite, name_main_test, name_mesh_short,
                   sprintf("%0*d", 4, n_nodes),
                   sprintf("%0*d", 4, n_locs),
                   sprintf("%0*d", 4, n_stat_units),
                   NSR_X_lc, NSR_Y, n_comp,
                   sep = "_"
                 )
                 json_list <- list(
                   test = list(
                     name_test = name_test
                   ),
                   mesh = list(
                     name_mesh = name_mesh
                   ),
                   dimensions = list(
                     n_nodes = n_nodes,
                     n_locs = n_locs,
                     n_stat_units = n_stat_units,
                     n_comp = n_comp,
                     n_reps = n_reps
                   ),
                   data = list(
                     locs_eq_nodes = locs_eq_nodes
                   ),
                   noise = list(
                     NSR_X_lc_vect = NSR_X_lc,
                     NSR_Y = NSR_Y
                   ),
                   COMPETITORS = FALSE
                 )
                 write_json(json_list, path = paste(path_queue, name_test, ".json", sep = ""))
               }
             }
           }
         },
         {
           stop(paste("The test", name_main_test, "does not exist"))
         }
  )
}
